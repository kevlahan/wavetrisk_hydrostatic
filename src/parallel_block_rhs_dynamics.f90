module parallel_block_rhs_dynamics_mod
  ! Direct connected RHS phases. Topology compilers and communication adapters
  ! bind addresses before execution; there is no expression interpreter here.
  use kind_mod, only : dp
  use parallel_block_rhs_mod, only : Rhs_Primitive_Stencil,rhs_column_start,rhs_pressure_patch, &
       rhs_primitives,rhs_pentagons,rhs_scalar_restrict,rhs_scalar_diffusion,rhs_velocity_diffusion
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  use parallel_block_rhs_diffusion_mod
  use parallel_block_mass_mod, only : execute_mass_restriction,execute_mass_divergence
  use parallel_block_velocity_mod, only : execute_velocity_sources,execute_velocity_gradients
  implicit none
  private
  public :: Rhs_Dynamics_Workspace,initialize_rhs_dynamics,rhs_dynamics_start,rhs_dynamics_layer
  public :: rhs_dynamics_diffusion,rhs_dynamics_primitives,rhs_dynamics_finish_level,rhs_dynamics_gradient
  public :: RHS_PRIMITIVE_VORT,RHS_FINAL_PV,RHS_SCALAR_FLUX,RHS_SCALAR_TENDENCY,RHS_VELOCITY_SOURCE
  public :: RHS_DIRECT_SCALAR_FLUX
  integer, parameter :: RHS_DIRECT_SCALAR_FLUX=20
  public :: RHS_BERNOULLI,RHS_EXNER,RHS_INJECTED_BERNOULLI,RHS_INJECTED_EXNER
  integer, parameter :: RHS_PRIMITIVE_VORT=11,RHS_FINAL_PV=12,RHS_SCALAR_FLUX=13,RHS_SCALAR_TENDENCY=14, &
       RHS_VELOCITY_SOURCE=15,RHS_BERNOULLI=16,RHS_EXNER=17,RHS_INJECTED_BERNOULLI=18,RHS_INJECTED_EXNER=19
  type :: Rhs_Dynamics_Workspace
     real(dp), allocatable :: lower(:),phi(:),pressure(:),midpoint(:),exner(:,:),thermodynamic(:,:)
     real(dp), allocatable :: flux(:,:),pv(:,:),ke(:),bernoulli(:),vorticity(:,:)
     real(dp), allocatable :: diffusion_flux(:,:),laplacian(:,:),divu(:,:),lap_divu(:,:),rotation(:,:)
     real(dp), allocatable :: source(:,:),physics(:,:),scalar_rhs(:,:),velocity_rhs(:,:)
     ! Physics values are generated immediately before their fused producer;
     ! SW values are indexed by event, never by an ambiguous destination alias.
     real(dp), allocatable :: scalar_physics(:,:),sw_physics(:,:)
  end type
  abstract interface
     subroutine Rhs_Phase_Halo(phase,level,values,context)
       import dp
       integer, intent(in) :: phase,level
       real(dp), intent(inout) :: values(:,:)
       class(*), intent(inout) :: context
     end subroutine
  end interface
contains
  subroutine initialize_rhs_dynamics(plan,nlayer,w)
    type(Rhs_Domain_Operators), intent(in) :: plan
    integer, intent(in) :: nlayer
    type(Rhs_Dynamics_Workspace), intent(out) :: w
    integer :: n,m,l
    if(.not.allocated(plan%dynamics)) error stop 'RHS dynamics program is absent'
    n=plan%nnode;m=0
    do l=lbound(plan%diffusion,1),ubound(plan%diffusion,1)
       m=max(m,size(plan%diffusion(l)%producer))
    end do
    allocate(w%lower(n),w%phi(n),w%pressure(n),w%midpoint(n),w%exner(n,nlayer),w%thermodynamic(n,2))
    allocate(w%flux(3*n,2),w%pv(3,n),w%ke(n),w%bernoulli(n),w%vorticity(2*n,1))
    allocate(w%diffusion_flux(3*n,2),w%laplacian(n,2),w%divu(n,1),w%lap_divu(n,1),w%rotation(3*n,1))
    allocate(w%source(3,n),w%physics(3,n),w%scalar_rhs(n,2),w%velocity_rhs(3,n))
    allocate(w%scalar_physics(3,m),w%sw_physics(3,m))
    w%lower=0;w%phi=0;w%pressure=0;w%midpoint=0;w%exner=0;w%thermodynamic=0
    w%flux=0;w%pv=0;w%ke=0;w%bernoulli=0;w%vorticity=0
    w%diffusion_flux=0;w%laplacian=0;w%divu=0;w%lap_divu=0;w%rotation=0
    w%source=0;w%physics=0;w%scalar_rhs=0;w%velocity_rhs=0
    w%scalar_physics=0;w%sw_physics=0
  end subroutine

  subroutine rhs_dynamics_start(plan,mass,mean_mass,surface,gravity,ptop,w)
    type(Rhs_Domain_Operators), intent(in) :: plan
    real(dp), intent(in) :: mass(:,:),mean_mass(:,:),surface(:),gravity,ptop
    type(Rhs_Dynamics_Workspace), intent(inout) :: w
    call rhs_column_start(plan%surface_visits,mass,mean_mass,surface,gravity,ptop,w%lower,w%phi)
  end subroutine

  subroutine rhs_dynamics_layer(plan,k,scalar,mean_scalar,velocity,order,nu,nu_pentagon, &
       gravity,cp,p0,kappa,w,halo,context)
    type(Rhs_Domain_Operators), intent(in) :: plan
    integer, intent(in) :: k,order(3)
    real(dp), intent(in) :: scalar(:,:),mean_scalar(:,:),velocity(:,:),nu(4),nu_pentagon,gravity,cp,p0,kappa
    type(Rhs_Dynamics_Workspace), intent(inout), target :: w
    procedure(Rhs_Phase_Halo) :: halo
    class(*), intent(inout) :: context
    integer :: l
    call rhs_dynamics_diffusion(plan,scalar,velocity,order,w,halo,context)
    w%source=0;w%scalar_rhs=0;w%velocity_rhs=0
    do l=ubound(plan%diffusion,1),lbound(plan%diffusion,1),-1
       call rhs_dynamics_primitives(plan,l,k,scalar,mean_scalar,velocity(:,1),order(1),nu(:2), &
            gravity,cp,p0,kappa,w,halo,context)
       call rhs_dynamics_finish_level(plan,l,k,order,nu,nu_pentagon,w,halo,context)
    end do
    call rhs_dynamics_gradient(plan,k,scalar,mean_scalar,w,halo,context)
  end subroutine

  subroutine rhs_dynamics_diffusion(plan,scalar,velocity,order,w,halo,context)
    type(Rhs_Domain_Operators), intent(in) :: plan
    real(dp), intent(in) :: scalar(:,:),velocity(:,:)
    integer, intent(in) :: order(3)
    type(Rhs_Dynamics_Workspace), intent(inout) :: w
    procedure(Rhs_Phase_Halo) :: halo
    class(*), intent(inout) :: context
    if(any(order<0).or.any(order>2)) error stop 'RHS unsupported diffusion order'
    if(size(scalar,2)/=2.or.size(velocity,2)/=1) error stop 'RHS layer variable extent'
    if(order(2)/=0) call rhs_divergence_sweep(plan%diffusion,velocity,w%diffusion_flux(:,:1),w%divu,halo,context)
    if(order(1)==2) call rhs_laplacian_sweep(plan%diffusion,scalar,w%diffusion_flux,w%laplacian,halo,context)
    if(order(2)==2) call rhs_laplacian_sweep(plan%diffusion,w%divu,w%diffusion_flux(:,:1),w%lap_divu,halo,context)
  end subroutine

  subroutine rhs_dynamics_primitives(plan,l,k,scalar,mean_scalar,velocity,order,nu,gravity,cp,p0,kappa,w,halo,context)
    type(Rhs_Domain_Operators), intent(in) :: plan
    integer, intent(in) :: l,k,order
    real(dp), intent(in) :: scalar(:,:),mean_scalar(:,:),velocity(:),nu(2),gravity,cp,p0,kappa
    type(Rhs_Dynamics_Workspace), intent(inout), target :: w
    procedure(Rhs_Phase_Halo) :: halo
    class(*), intent(inout) :: context
    real(dp), pointer :: mass_flux(:,:),pvview(:,:)
    integer :: j,a,b,c,d
    mass_flux(1:3,1:plan%nnode)=>w%flux(:,1)
    do j=1,size(plan%pressure(l)%patch)
       associate(p=>plan%pressure(l)%patch(j))
       a=p%first_pressure;b=p%last_pressure;c=p%first_producer;d=p%last_producer
       call rhs_pressure_patch(plan%pressure(l)%node(a:b),scalar(:,1),scalar(:,2), &
            mean_scalar(:,1),mean_scalar(:,2),gravity,cp,p0,kappa,w%lower,w%phi,w%pressure,w%exner(:,k),w%midpoint)
       call scalar_physics(plan%diffusion(l)%producer(c:d),scalar(:,1),w%laplacian(:,1),nu(1),order, &
            w%scalar_physics(:,:d-c+1),w%sw_physics(:,:d-c+1))
       call rhs_primitives(plan%diffusion(l)%producer(c:d),scalar(:,1),mean_scalar(:,1),velocity,w%midpoint, &
            w%physics,mass_flux,w%pv,w%ke,w%bernoulli,w%vorticity(:,1),w%sw_physics(:,:d-c+1), &
            event_physics=w%scalar_physics(:,:d-c+1))
       call scalar_physics(plan%diffusion(l)%producer(c:d),scalar(:,2),w%laplacian(:,2),nu(2),order, &
            w%scalar_physics(:,:d-c+1),w%sw_physics(:,:d-c+1))
       call temperature_flux(plan%diffusion(l)%producer(c:d),scalar(:,2),mean_scalar(:,2),velocity, &
            w%scalar_physics(:,:d-c+1),w%sw_physics(:,:d-c+1),w%flux(:,2))
       end associate
    end do
    call halo(RHS_PRIMITIVE_VORT,l,w%vorticity,context)
    call rhs_pentagons(plan%rotation(l)%pentagon,scalar(:,1),mean_scalar(:,1),velocity,w%pv,w%vorticity(:,1))
    ! Separate PV and triangle readiness: pentagons write different field sets.
    pvview(1:3*plan%nnode,1:1)=>w%pv
    call halo(RHS_FINAL_PV,l,pvview,context)
    call halo(RHS_ROTATION_TRIANGLES,l,w%vorticity,context)
    ! Both fields have the same canonical producer/reader identities and are
    ! complete here. Pack them into one persistent two-column transaction.
    w%thermodynamic(:,1)=w%bernoulli;w%thermodynamic(:,2)=w%exner(:,k)
    call halo(RHS_BERNOULLI,l,w%thermodynamic,context)
    w%bernoulli=w%thermodynamic(:,1);w%exner(:,k)=w%thermodynamic(:,2)
    call rhs_scalar_restrict(plan%dynamics(l)%injection(1,:),plan%dynamics(l)%injection(2,:),w%bernoulli,w%exner(:,k))
    w%thermodynamic(:,1)=w%bernoulli;w%thermodynamic(:,2)=w%exner(:,k)
    call halo(RHS_INJECTED_BERNOULLI,l,w%thermodynamic,context)
    w%bernoulli=w%thermodynamic(:,1);w%exner(:,k)=w%thermodynamic(:,2)
  end subroutine

  subroutine scalar_physics(plan,scalar,laplacian,nu,order,physics,sw)
    type(Rhs_Primitive_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: scalar(:),laplacian(:),nu
    integer, intent(in) :: order
    real(dp), intent(inout) :: physics(:,:)
    real(dp), intent(out) :: sw(:,:)
    integer :: s,n(4)
    real(dp) :: x(4)
    do s=1,size(plan)
       associate(p=>plan(s))
       n=p%mass([1,3,6,2]);x=0
       if(order==1) x=scalar(n)
       if(order==2) x=laplacian(n)
       physics(:,s)=rhs_scalar_diffusion(x,p%geometry([1,4,8]),p%geometry(11:13),nu,order)
       sw(:,s)=0
       if(p%southwest) then
          n=[p%mass(1),p%mass(5),p%sw_mass,p%mass(4)]
          if(order==1) x=scalar(n)
          if(order==2) x=laplacian(n)
          sw(:,s)=rhs_scalar_diffusion(x,-p%geometry([3,6,10]),p%geometry(14:16),nu,order)
       end if
       end associate
    end do
  end subroutine

  subroutine temperature_flux(plan,scalar,mean,velocity,physics,sw,flux)
    type(Rhs_Primitive_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: scalar(:),mean(:),velocity(:),physics(:,:),sw(:,:)
    real(dp), intent(inout) :: flux(:)
    real(dp) :: x(4),v(3)
    integer :: s,n(4),t,e
    do s=1,size(plan)
       associate(p=>plan(s))
       n=p%mass([1,3,6,2]);x=scalar(n)+mean(n)
       v=velocity(p%velocity([1,4,8]))*p%geometry(11:13)
       t=3*(p%target-1)
       flux(t+1:t+3)=v*(0.5_dp*(x(1)+x(2:)))+physics(:,s)
       if(p%southwest) then
          n=[p%mass(1),p%mass(5),p%sw_mass,p%mass(4)];x=scalar(n)+mean(n)
          v=velocity(p%velocity([3,6,10]))*p%geometry(14:16)
          v=v*(0.5_dp*(x(1)+x(2:)))+sw(:,s)
          do e=1,3
             t=3*(p%sw_edge(2,e)-1)+p%sw_edge(1,e);flux(t)=v(e)
          end do
       end if
       end associate
    end do
  end subroutine

  subroutine rhs_dynamics_finish_level(plan,l,k,order,nu,nu_pentagon,w,halo,context)
    type(Rhs_Domain_Operators), intent(in) :: plan
    integer, intent(in) :: l,k,order(3)
    real(dp), intent(in) :: nu(4),nu_pentagon
    type(Rhs_Dynamics_Workspace), intent(inout), target :: w
    procedure(Rhs_Phase_Halo) :: halo
    class(*), intent(inout) :: context
    real(dp), pointer :: mass_flux(:,:),sourceview(:,:)
    real(dp) :: div(4),vort(4),coefficient
    integer :: v,s,n
    call halo(RHS_DIRECT_SCALAR_FLUX,l,w%flux,context)
    do v=1,2
       call execute_mass_restriction(plan%diffusion(l)%restriction,w%flux(:,v),w%scalar_rhs(:,v),record_work=.false.)
    end do
    call halo(RHS_SCALAR_FLUX,l,w%flux,context)
    do v=1,2
       call execute_mass_divergence(plan%dynamics(l)%divergence,w%flux(:,v),w%scalar_rhs(:,v),record_work=.false.)
    end do
    if(l>lbound(plan%diffusion,1)) call halo(RHS_SCALAR_TENDENCY,l,w%scalar_rhs,context)
    if(order(3)==2) call rhs_rotational_laplacian(l,plan%rotation(l)%curl,plan%diffusion(l)%producer, &
         plan%rotation(l)%pentagon,w%vorticity,w%rotation,halo,context)
    do s=1,size(plan%dynamics(l)%physics)
       associate(p=>plan%dynamics(l)%physics(s))
       n=p%node(1);div=0;vort=0
       if(order(2)==1) div=w%divu(p%node,1)
       if(order(2)==2) div=w%lap_divu(p%node,1)
       if(order(3)/=0) vort=w%vorticity(p%vorticity,1)
       coefficient=nu(3)
       if(order(2)==1.and.p%pentagon) coefficient=nu_pentagon
       w%physics(:,n)=rhs_velocity_diffusion(div,vort,p%length,p%dual_length,coefficient,nu(4),order(2),order(3))
       end associate
    end do
    mass_flux(1:3,1:plan%nnode)=>w%flux(:,1)
    call execute_velocity_sources(plan%dynamics(l)%source,mass_flux,w%pv,w%physics,w%source)
    sourceview(1:3*plan%nnode,1:1)=>w%source
    call halo(RHS_VELOCITY_SOURCE,l,sourceview,context)
  end subroutine

  subroutine rhs_dynamics_gradient(plan,k,scalar,mean_scalar,w,halo,context)
    type(Rhs_Domain_Operators), intent(in) :: plan
    integer, intent(in) :: k
    real(dp), intent(in) :: scalar(:,:),mean_scalar(:,:)
    type(Rhs_Dynamics_Workspace), intent(inout) :: w
    procedure(Rhs_Phase_Halo) :: halo
    class(*), intent(inout) :: context
    ! Caller binds final B/Exner readers once, after all fine-to-coarse injection.
    call execute_velocity_gradients(plan%gradient,w%bernoulli,w%exner(:,k),mean_scalar(:,1)+scalar(:,1), &
         mean_scalar(:,2)+scalar(:,2),w%source,w%velocity_rhs)
  end subroutine
end module
