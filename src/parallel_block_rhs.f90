module parallel_block_rhs_mod
  ! Connected numerical components for the block-owner RHS cutover.
  ! Development only: not selected by production until routes and RK callers move.
  ! All addresses are one-based local workspace slots, including explicit halos.
  use kind_mod, only : dp
  use parallel_block_mass_mod, only : Mass_Level, execute_mass_restriction, execute_mass_divergence
  use parallel_block_velocity_mod, only : Velocity_Level_Program, Velocity_Gradient_Stencil, &
       execute_velocity_sources, execute_velocity_gradients, velocity_pressure_layer
  implicit none
  private
  public :: Rhs_Pentagon_Stencil, rhs_pentagons
  public :: Rhs_Primitive_Stencil, rhs_columns, rhs_primitives
  public :: rhs_column_start,rhs_pressure_patch
  public :: rhs_scalar_restrict, rhs_scalar_diffusion, rhs_velocity_diffusion
  public :: rhs_mass_restrict, rhs_mass_divergence, rhs_velocity_source, rhs_velocity_gradient, rhs_rk

  ! Same C,N,E,S,W,NE mass and velocity mapping as Stage276/277.
  ! The compiler must preserve step1 traversal order, including repeated vort writes.
  type :: Rhs_Primitive_Stencil
     integer :: target=0, mass(6)=0, velocity(10)=0, vorticity(4)=0
     real(dp) :: geometry(37)=0.0_dp
     ! comp_SW follows this node's comput immediately, before the next node.
     ! Edge destinations are (component,node), already resolved by the compiler.
     logical :: southwest=.false.
     integer :: sw_mass=0,sw_velocity(2)=0,sw_edge(2,3)=0,sw_vorticity(2)=0
     ! RT/UP SW lengths; six mass-area parts; two Coriolis and triangle areas.
     real(dp) :: sw_geometry(12)=0.0_dp
  end type
  ! Uses established shared_mod corner values: IJPLUS5, IPLUSJMINUS6,
  ! IJMINUS7, IMINUSJPLUS8. Nodes are C,N,E,S,W,NE,SW.
  type :: Rhs_Pentagon_Stencil
     integer :: corner=0,mass(7)=0,velocity(3,7)=0,vorticity(2,7)=0,pv_node(7)=0
     real(dp) :: length(3,7)=0.0_dp,triangle_area(2,7)=0.0_dp,coriolis(2,7)=0.0_dp,parts(6,7)=0.0_dp
  end type
  integer, parameter :: IJPLUS=5,IPLUSJMINUS=6,IJMINUS=7,IMINUSJPLUS=8
  integer, parameter :: dual(6)=[1,4,8,3,6,10]
contains
  subroutine rhs_column_start(visits,mass,mean_mass,surface,gravity,ptop,lower,geopotential)
    ! Surface aliases may be visited more than once by the accepted patch
    ! traversal. Repeated additions are part of its floating-point contract.
    integer, intent(in) :: visits(:)
    real(dp), intent(in) :: mass(:,:),mean_mass(:,:),surface(:),gravity,ptop
    real(dp), intent(out) :: lower(:)
    real(dp), intent(inout) :: geopotential(:)
    integer :: n,k,j
    if(size(visits)/=size(mass,1).or.any(shape(mass)/=shape(mean_mass)).or. &
         size(lower)/=size(visits).or.size(surface)/=size(visits).or.size(geopotential)/=size(visits)) &
         error stop 'RHS column start extent'
    if(any(visits<0)) error stop 'RHS column negative visits'
    lower=0.0_dp
    do k=1,size(mass,2)
       do n=1,size(visits)
          do j=1,visits(n)
             lower(n)=lower(n)+(mean_mass(n,k)+mass(n,k))
          end do
       end do
    end do
    lower=gravity*lower+ptop
    where(visits>0) geopotential=surface
  end subroutine

  subroutine rhs_pressure_patch(nodes,mass,temperature,mean_mass,mean_temperature, &
       gravity,cp,p0,kappa,lower,geopotential,pressure,exner,midpoint)
    ! Call immediately before this patch's primitives, preserving repeated
    ! physical alias updates within a layer and the fine-to-coarse level order.
    integer, intent(in) :: nodes(:)
    real(dp), intent(in) :: mass(:),temperature(:),mean_mass(:),mean_temperature(:),gravity,cp,p0,kappa
    real(dp), intent(inout) :: lower(:),geopotential(:),pressure(:),exner(:),midpoint(:)
    integer :: j,n
    real(dp) :: upper,phi
    do j=1,size(nodes)
       n=nodes(j)
       call velocity_pressure_layer(mean_mass(n)+mass(n),mean_temperature(n)+temperature(n), &
            gravity,cp,p0,kappa,lower(n),geopotential(n),pressure(n),upper,exner(n),phi)
       midpoint(n)=0.5_dp*(phi+geopotential(n))
       lower(n)=upper;geopotential(n)=phi
    end do
  end subroutine

  subroutine rhs_columns(mass,temperature,mean_mass,mean_temperature,surface_geopotential, &
       gravity,cp,p0,kappa,ptop,pressure,exner,geopotential_midpoint)
    ! Complete physical columns on the executing owner. State and means stay
    ! separate: pressure uses mean+state, while step1 uses state+mean.
    real(dp), intent(in) :: mass(:,:),temperature(:,:),mean_mass(:,:),mean_temperature(:,:)
    real(dp), intent(in) :: surface_geopotential(:),gravity,cp,p0,kappa,ptop
    real(dp), intent(out) :: pressure(:,:),exner(:,:),geopotential_midpoint(:,:)
    integer :: n,k,nnode,nlayer
    real(dp) :: lower_pressure,upper_pressure,lower_geopotential,upper_geopotential,rho,theta
    nnode=size(mass,1);nlayer=size(mass,2)
    if (any(shape(temperature)/=shape(mass)).or.any(shape(mean_mass)/=shape(mass)).or. &
         any(shape(mean_temperature)/=shape(mass)).or.size(surface_geopotential)/=nnode.or. &
         any(shape(pressure)/=shape(mass)).or.any(shape(exner)/=shape(mass)).or. &
         any(shape(geopotential_midpoint)/=shape(mass))) error stop "block RHS column extents differ"
    do n=1,nnode
       lower_pressure=0.0_dp
       do k=1,nlayer
          rho=mean_mass(n,k)+mass(n,k)
          lower_pressure=lower_pressure+rho
       end do
       lower_pressure=gravity*lower_pressure+ptop
       lower_geopotential=surface_geopotential(n)
       do k=1,nlayer
          rho=mean_mass(n,k)+mass(n,k)
          theta=mean_temperature(n,k)+temperature(n,k)
          call velocity_pressure_layer(rho,theta,gravity,cp,p0,kappa,lower_pressure,lower_geopotential, &
               pressure(n,k),upper_pressure,exner(n,k),upper_geopotential)
          geopotential_midpoint(n,k)=0.5_dp*(upper_geopotential+lower_geopotential)
          lower_pressure=upper_pressure
          lower_geopotential=upper_geopotential
       end do
    end do
  end subroutine

  subroutine rhs_primitives(plan,mass,mean_mass,velocity,geopotential_midpoint,physics_flux, &
       flux,pv,ke,bernoulli,vorticity,sw_physics,event_physics)
    ! Retain direct fusion. Physics is a named phase input, not a Domain flux.
    ! No publication of the 19-value validation packet is needed in production.
    type(Rhs_Primitive_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: mass(:),mean_mass(:),velocity(:),geopotential_midpoint(:),physics_flux(:,:)
    real(dp), intent(inout) :: flux(:,:),pv(:,:),ke(:),bernoulli(:),vorticity(:)
    real(dp), optional, intent(in) :: sw_physics(:,:),event_physics(:,:)
    integer :: s,t,e,c,n
    real(dp) :: x(16),y(19),sw(8),sw_mass,sw_velocity(2)
    if (size(mass)/=size(mean_mass).or.size(flux,1)/=3.or.any(shape(pv)/=shape(flux)).or. &
         any(shape(physics_flux)/=shape(flux)).or.size(ke)/=size(flux,2).or. &
         size(bernoulli)/=size(ke).or.size(geopotential_midpoint)/=size(ke)) &
         error stop "block RHS primitive extents differ"
    if(present(event_physics)) then
       if(any(shape(event_physics)/=[3,size(plan)])) error stop 'RHS event physics extent'
    end if
    ! Reject an incomplete SW transaction before any output writes.
    if (any(plan%southwest)) then
       if (.not.present(sw_physics)) error stop "block RHS SW physics missing"
       if (any(shape(sw_physics)/=[3,size(plan)])) error stop "block RHS SW physics extents differ"
       do s=1,size(plan)
          associate(p=>plan(s))
          if (.not.p%southwest) cycle
          if (p%sw_mass<1.or.p%sw_mass>size(mass).or.any(p%sw_velocity<1).or. &
               any(p%sw_velocity>size(velocity)).or.any(p%sw_edge(1,:)<1).or.any(p%sw_edge(1,:)>3).or. &
               any(p%sw_edge(2,:)<1).or.any(p%sw_edge(2,:)>size(flux,2)).or. &
               any(p%sw_vorticity<1).or.any(p%sw_vorticity>size(vorticity))) &
               error stop "block RHS SW stencil invalid"
          end associate
       end do
    end if
    do s=1,size(plan)
       associate(p=>plan(s))
       t=p%target
       if (t<1.or.t>size(ke).or.any(p%mass<1).or.any(p%mass>size(mass)).or. &
            any(p%velocity<1).or.any(p%velocity>size(velocity)).or. &
            any(p%vorticity<1).or.any(p%vorticity>size(vorticity))) &
            error stop "block RHS primitive stencil invalid"
       x(1:6)=mass(p%mass)+mean_mass(p%mass)
       x(7:16)=velocity(p%velocity)
       call fused_primitive(x,p%geometry,y)
       if(present(event_physics)) then
          flux(:,t)=y(1:3)+event_physics(:,s)
       else
          flux(:,t)=y(1:3)+physics_flux(:,t)
       end if
       pv(:,t)=y(4:6)
       ke(t)=y(7)
       bernoulli(t)=ke(t)+geopotential_midpoint(t)
       ! Explicit sequential assignments preserve final shared-vorticity writes.
       vorticity(p%vorticity(1))=y(8)
       vorticity(p%vorticity(2))=y(9)
       vorticity(p%vorticity(3))=y(10)
       vorticity(p%vorticity(4))=y(11)
       if (p%southwest) then
          sw_mass=mass(p%sw_mass)+mean_mass(p%sw_mass)
          sw_velocity=velocity(p%sw_velocity)
          call fused_southwest(x,y,p%geometry,p%sw_geometry,sw_mass,sw_velocity,sw)
          ! Same qe, triangle, then flux order as comp_SW. Aliases retain
          ! sequential final-writer semantics; no regrouping into field passes.
          do e=1,3
             c=p%sw_edge(1,e);n=p%sw_edge(2,e)
             pv(c,n)=sw(3+e)
          end do
          vorticity(p%sw_vorticity(1))=sw(7)
          vorticity(p%sw_vorticity(2))=sw(8)
          do e=1,3
             c=p%sw_edge(1,e);n=p%sw_edge(2,e)
             flux(c,n)=sw(e)+sw_physics(e,s)
          end do
       end if
       end associate
    end do
  end subroutine

  pure subroutine fused_primitive(x,g,y)
    ! Verified Stage277 arithmetic, without its geometry-owner cache.
    real(dp), intent(in) :: x(16),g(37)
    real(dp), intent(out) :: y(19)
    real(dp) :: u(10),v(6),c(4),m(4),pv(4)
    u=x(7:16)*g(1:10)
    v=x(6+dual)*g(11:16)
    c(1)=u(1)+u(9)+u(4);c(2)=-(u(4)+u(8)+u(2))
    c(3)=u(3)+u(8)+u(7);c(4)=-(u(1)+u(5)+u(10))
    m(1)=x(1)*g(17)+x(3)*g(18)+x(6)*g(19)
    m(2)=x(1)*g(20)+x(6)*g(21)+x(2)*g(22)
    m(3)=x(5)*g(23)+x(1)*g(24)+x(2)*g(25)
    m(4)=x(4)*g(26)+x(3)*g(27)+x(1)*g(28)
    pv=(g(29:32)+c)/m
    y(1)=v(1)*(0.5_dp*(x(1)+x(3)))
    y(2)=v(2)*(0.5_dp*(x(1)+x(6)))
    y(3)=v(3)*(0.5_dp*(x(1)+x(2)))
    y(4)=0.5_dp*(pv(4)+pv(1));y(5)=0.5_dp*(pv(2)+pv(1));y(6)=0.5_dp*(pv(2)+pv(3))
    y(7)=(u(8)*v(3)+u(4)*v(2)+u(1)*v(1)+u(10)*v(6)+u(6)*v(5)+u(3)*v(4))*g(33)/4
    y(8:11)=c/g(34:37)
    y(12:13)=pv(3:4)
    y(14:19)=v
  end subroutine

  pure subroutine fused_southwest(x,direct,g,geometry,mass,velocity,value)
    ! Direct comp_SW arithmetic; direct retains this comput's corner PV and
    ! dual velocities. Recomputing corner PV from stored qe would be incorrect.
    real(dp), intent(in) :: x(16),direct(19),g(37),geometry(12),mass,velocity(2)
    real(dp), intent(out) :: value(8)
    real(dp) :: u_rt_w,u_dg_sw,u_up_s,u_rt_sw,u_up_sw,circ(2),pv_mass(2),pv(2)
    u_rt_w=x(9)*g(3)
    u_dg_sw=x(12)*g(6)
    u_up_s=x(16)*g(10)
    u_rt_sw=velocity(1)*geometry(1)
    u_up_sw=velocity(2)*geometry(2)
    circ(1)=u_rt_sw+u_up_s+u_dg_sw
    circ(2)=-(u_rt_w+u_dg_sw+u_up_sw)
    pv_mass(1)=mass*geometry(3)+x(4)*geometry(4)+x(1)*geometry(5)
    pv_mass(2)=mass*geometry(6)+x(1)*geometry(7)+x(5)*geometry(8)
    pv=(geometry(9:10)+circ)/pv_mass
    value(1)=direct(17)*(0.5_dp*(x(1)+x(5)))
    value(2)=direct(18)*(0.5_dp*(x(1)+mass))
    value(3)=direct(19)*(0.5_dp*(x(1)+x(4)))
    value(4)=0.5_dp*(direct(12)+pv(2))
    value(5)=0.5_dp*(pv(1)+pv(2))
    value(6)=0.5_dp*(pv(1)+direct(13))
    value(7:8)=circ/geometry(11:12)
  end subroutine

  subroutine rhs_pentagons(plan,mass,mean_mass,velocity,pv,vorticity)
    ! Run after direct/SW producers and the required final-writer halo.
    ! Mapped writes preserve post_step1 order, including actual slot aliases.
    type(Rhs_Pentagon_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: mass(:),mean_mass(:),velocity(:)
    real(dp), intent(inout) :: pv(:,:),vorticity(:)
    integer :: s
    real(dp) :: rho(7),pv_mass(3)
    real(dp) :: circ_LORT,circ_UPLT,circ_S_UPLT,circ_SW_LORT,circ_SW_UPLT,circ_W_LORT
    real(dp) :: pv_LORT,pv_UPLT,pv_S_UPLT,pv_SW_LORT,pv_SW_UPLT,pv_W_LORT
    real(dp) :: u_prim_RT,u_prim_RT_N,u_prim_RT_SW,u_prim_RT_W
    real(dp) :: u_prim_DG_SW,u_prim_UP,u_prim_UP_S,u_prim_UP_SW
    if (size(mass)/=size(mean_mass).or.size(pv,1)/=3) error stop "block RHS pentagon extents differ"
    do s=1,size(plan)
       associate(p=>plan(s))
       if (p%corner<IJPLUS.or.p%corner>IMINUSJPLUS.or.any(p%mass<1).or.any(p%mass>size(mass)).or. &
            any(p%velocity<1).or.any(p%velocity>size(velocity)).or. &
            any(p%vorticity<1).or.any(p%vorticity>size(vorticity)).or. &
            any(p%pv_node<1).or.any(p%pv_node>size(pv,2))) error stop "block RHS pentagon stencil invalid"
       end associate
    end do
    do s=1,size(plan)
       associate(p=>plan(s))
       rho=mass(p%mass)+mean_mass(p%mass)
    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (p%corner == IJMINUS) then


       u_prim_RT_W  = velocity(p%velocity(1,5)) * p%length(1,5)
       u_prim_RT_SW = velocity(p%velocity(1,7)) * p%length(1,7) 
       u_prim_UP_S  = velocity(p%velocity(3,4)) * p%length(3,4)

       circ_SW_LORT = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW
       circ_W_LORT  = vorticity(p%vorticity(1,5)) * p%triangle_area(1,5)
       circ_S_UPLT  = vorticity(p%vorticity(2,4)) * p%triangle_area(2,4)

       vorticity(p%vorticity(1,7)) = circ_SW_LORT / p%triangle_area(1,7)
       vorticity(p%vorticity(2,7)) = vorticity(p%vorticity(1,7))

       pv_mass(1) = &
            rho(5)  * p%parts(6,5) + &
            rho(1)     * sum(p%parts(4:5,1)) + &
            rho(4) * p%parts(3,4)
      
       pv_mass(2) = &
            rho(5)  * p%parts(1,5) + &
            rho(1)     * p%parts(3,1) + &
            rho(2) * p%parts(5,2)
       
       pv_mass(3) = &
            rho(4) * p%parts(2,4) + &
            rho(3)  * p%parts(4,3) + &
            rho(1)     * p%parts(6,1)

       pv_W_LORT  = (p%coriolis(1,5) + circ_W_LORT ) / pv_mass(1)
       pv_SW_LORT = (p%coriolis(1,7) + circ_SW_LORT) / pv_mass(2)
       pv_S_UPLT  = (p%coriolis(2,4) + circ_S_UPLT ) / pv_mass(3)

       pv_SW_UPLT = pv_SW_LORT

       pv(1,p%pv_node(5)) = (0.5_dp*(pv_W_LORT+ pv_SW_UPLT))
       pv(3,p%pv_node(4)) = (0.5_dp*(pv_S_UPLT+ pv_SW_LORT))
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
    if (p%corner == IPLUSJMINUS) then 


       u_prim_RT_SW = velocity(p%velocity(1,7)) * p%length(1,7)
       u_prim_DG_SW = velocity(p%velocity(2,7)) * p%length(2,7)
       u_prim_RT    = velocity(p%velocity(1,1)) * p%length(1,1)

       circ_SW_LORT = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW
       circ_LORT    = vorticity(p%vorticity(1,1)) * p%triangle_area(1,1)
       circ_SW_UPLT = vorticity(p%vorticity(2,7)) * p%triangle_area(2,7)
       
       vorticity(p%vorticity(1,7)) = circ_SW_LORT / p%triangle_area(1,7)
       vorticity(p%vorticity(2,4)) = vorticity(p%vorticity(1,7))

       pv_mass(1) = &
            rho(7) * p%parts(1,7) + &
            rho(4)     * p%parts(3,4) + &
            rho(1)         * sum(p%parts(5:6,1))

       pv_mass(2) = &
            rho(1)         * p%parts(1,1) + &
            rho(3)      * p%parts(3,3) + &
            rho(6) * p%parts(5,6)

       pv_mass(3) = &
            rho(7) * p%parts(2,7) + &
            rho(1)         * p%parts(4,1) + &
            rho(5)      * p%parts(6,5)

       pv_SW_LORT = (p%coriolis(1,7) + circ_SW_LORT) / pv_mass(1)
       pv_LORT    = (p%coriolis(1,1) + circ_LORT   ) / pv_mass(2)
       pv_SW_UPLT = (p%coriolis(2,7) + circ_SW_UPLT) / pv_mass(3)

       pv_S_UPLT = pv_SW_LORT

       pv(1,p%pv_node(1)) = (0.5_dp*(pv_S_UPLT+  pv_LORT))
       pv(2,p%pv_node(7)) = (0.5_dp*(pv_SW_LORT+ pv_SW_UPLT))
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
    if (p%corner == IMINUSJPLUS) then


       u_prim_UP    = velocity(p%velocity(3,1)) * p%length(3,1)
       u_prim_DG_SW = velocity(p%velocity(2,7)) * p%length(2,7)
       u_prim_UP_SW = velocity(p%velocity(3,7)) * p%length(3,7)

       circ_SW_UPLT = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW
       circ_UPLT    = vorticity(p%vorticity(2,1)) * p%triangle_area(2,1)
       circ_SW_LORT = vorticity(p%vorticity(1,7)) * p%triangle_area(1,7)

       vorticity(p%vorticity(2,7)) = circ_SW_UPLT / p%triangle_area(2,7)  
       vorticity(p%vorticity(1,5)) = vorticity(p%vorticity(2,7))

       pv_mass(1) = &
            rho(7) * p%parts(2,7) + &
            rho(1)         * sum(p%parts(3:4,1)) + &
            rho(5)      * p%parts(6,5)

       pv_mass(2) = &
            rho(1)         * p%parts(2,1) + &
            rho(6) * p%parts(4,6) + &
            rho(2)     * p%parts(6,2)

       pv_mass(3) = &
            rho(7) * p%parts(1,7) + &
            rho(4)     * p%parts(3,4) + &
            rho(1)         * p%parts(5,1)
       
       pv_SW_UPLT = (p%coriolis(2,7) + circ_SW_UPLT) / pv_mass(1)
       pv_UPLT    = (p%coriolis(2,1) + circ_UPLT   ) / pv_mass(2)
       pv_SW_LORT = (p%coriolis(1,7) + circ_SW_LORT) / pv_mass(3)           

       pv_W_LORT = pv_SW_UPLT

       pv(3,p%pv_node(1)) = (0.5_dp*(pv_W_LORT+  pv_UPLT))
       pv(2,p%pv_node(7)) = (0.5_dp*(pv_SW_LORT+ pv_SW_UPLT))
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*id+DG+1) = 0 in this case
    if (p%corner == IJPLUS) then 


       u_prim_RT   = velocity(p%velocity(1,1)) * p%length(1,1)
       u_prim_RT_N = velocity(p%velocity(1,2)) * p%length(1,2)
       u_prim_UP   = velocity(p%velocity(3,1)) * p%length(3,1)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP
       circ_W_LORT = vorticity(p%vorticity(1,5)) * p%triangle_area(1,5)
       circ_S_UPLT = vorticity(p%vorticity(2,4)) * p%triangle_area(2,4)

       vorticity(p%vorticity(1,1)) = circ_LORT / p%triangle_area(1,1)
       vorticity(p%vorticity(2,1)) = vorticity(p%vorticity(1,1))

       pv_mass(1) = &
            rho(3)  * p%parts(3,3) + &
            rho(1)     * sum(p%parts(1:2,1)) + &
            rho(2) * p%parts(6,2)

       pv_mass(2) = &
            rho(5)  * p%parts(1,5) + &
            rho(1)     * p%parts(3,1) + &
            rho(2) * p%parts(5,2)

       pv_mass(3) = &
            rho(4) * p%parts(2,4) + &
            rho(3)  * p%parts(4,3) + &
            rho(1)     * p%parts(6,1)

       pv_LORT   = (p%coriolis(1,1) + circ_LORT  ) / pv_mass(1)                   
       pv_W_LORT = (p%coriolis(1,5) + circ_W_LORT) / pv_mass(2)
       pv_S_UPLT = (p%coriolis(2,4) + circ_S_UPLT) / pv_mass(3)

       pv_UPLT = pv_LORT

       pv(1,p%pv_node(1)) = (0.5_dp*(pv_LORT+ pv_S_UPLT))
       pv(3,p%pv_node(1)) = (0.5_dp*(pv_UPLT+ pv_W_LORT))
    end if
       end associate
    end do
  end subroutine rhs_pentagons

  subroutine rhs_scalar_restrict(target,child,bernoulli,exner)
    ! cpt_or_restr_scalar is an ordered injection, not an averaging reduction.
    ! Compiler emits only mask_n>=RESTRCT writes, preserving child traversal.
    integer, intent(in) :: target(:),child(:)
    real(dp), intent(inout) :: bernoulli(:),exner(:)
    integer :: s
    if (size(target)/=size(child).or.size(bernoulli)/=size(exner)) &
         error stop "block RHS scalar restriction extents differ"
    if (any(target<1).or.any(target>size(bernoulli)).or.any(child<1).or.any(child>size(bernoulli))) &
         error stop "block RHS scalar restriction address invalid"
    ! Legacy processes the entire Bernoulli sweep before the Exner sweep.
    do s=1,size(target)
       bernoulli(target(s))=bernoulli(child(s))
    end do
    do s=1,size(target)
       exner(target(s))=exner(child(s))
    end do
  end subroutine

  pure function rhs_scalar_diffusion(scalar,length,dual_length,nu,order) result(flux)
    ! climate physics_scalar_flux_case. For order 2 scalar is the completed
    ! Laplacian, never state. Signed length encodes the SW adapter orientation.
    real(dp), intent(in) :: scalar(4),length(3),dual_length(3),nu
    integer, intent(in) :: order
    real(dp) :: flux(3),grad(3)
    flux=0.0_dp
    if (order==0) return
    grad=[(scalar(2)-scalar(1))/length(1),(scalar(1)-scalar(3))/length(2), &
         (scalar(4)-scalar(1))/length(3)]
    flux=(-1)**order*nu*grad*dual_length
  end function

  pure function rhs_velocity_diffusion(divergence,vorticity,length,dual_length,nu_div,nu_rot, &
       order_div,order_rot) result(source)
    ! climate physics_velo_source_case. Caller supplies pentagon-adjusted nu_div
    ! and final post_step1/post_vort values; no shared C_visc mutation is needed.
    real(dp), intent(in) :: divergence(4),vorticity(4),length(3),dual_length(3),nu_div,nu_rot
    integer, intent(in) :: order_div,order_rot
    real(dp) :: source(3),grad(3),curl(3)
    source=0.0_dp
    if (order_div/=0) then
       grad=[(divergence(2)-divergence(1))/length(1),(divergence(1)-divergence(3))/length(2), &
            (divergence(4)-divergence(1))/length(3)]
       source=(-1)**(order_div-1)*nu_div*grad
    end if
    if (order_rot/=0) then
       ! Vorticity inputs are C-LORT, C-UPLT, W-LORT, S-UPLT.
       curl=[(vorticity(1)-vorticity(4))/dual_length(1), &
            (vorticity(1)-vorticity(2))/dual_length(2),(vorticity(3)-vorticity(2))/dual_length(3)]
       source=source-(-1)**(order_rot-1)*nu_rot*curl
    end if
  end function

  subroutine rhs_mass_restrict(level,flux,tendency)
    ! Requires completed fine-level node halo. Edge halo follows this call.
    type(Mass_Level), intent(in) :: level
    real(dp), intent(inout) :: flux(:)
    real(dp), intent(in) :: tendency(:)
    call execute_mass_restriction(level%restriction,flux,tendency)
  end subroutine

  subroutine rhs_mass_divergence(level,flux,tendency)
    ! Requires restricted edge halo. Node halo precedes the next coarser level.
    type(Mass_Level), intent(in) :: level
    real(dp), intent(in) :: flux(:)
    real(dp), intent(inout) :: tendency(:)
    call execute_mass_divergence(level%divergence,flux,tendency)
  end subroutine

  subroutine rhs_velocity_source(level,flux,pv,physics,source)
    ! Consume the same restricted flux as mass divergence. Compile ordered
    ! direct/restricted actions into local slots; exchange between levels only.
    type(Velocity_Level_Program), intent(in) :: level
    real(dp), intent(in) :: flux(:,:),pv(:,:),physics(:,:)
    real(dp), intent(inout) :: source(:,:)
    call execute_velocity_sources(level,flux,pv,physics,source)
  end subroutine

  subroutine rhs_velocity_gradient(plan,bernoulli,exner,mass,temperature,source,tendency)
    ! FINAL restricted B/Exner, never independently reconstructed gradients.
    type(Velocity_Gradient_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: bernoulli(:),exner(:),mass(:),temperature(:),source(:,:)
    real(dp), intent(inout) :: tendency(:,:)
    call execute_velocity_gradients(plan,bernoulli,exner,mass,temperature,source,tendency)
  end subroutine

  subroutine rhs_rk(base_mass,base_velocity,mass_rhs,velocity_rhs,h,stage_mass,stage_velocity)
    ! Both authoritative fields consume the connected RHS. Outputs must not
    ! alias inputs; final RK stage uses its saved resident base state.
    real(dp), intent(in) :: base_mass(:),base_velocity(:,:),mass_rhs(:),velocity_rhs(:,:),h
    real(dp), intent(out) :: stage_mass(:),stage_velocity(:,:)
    if (size(base_mass)/=size(mass_rhs).or.size(stage_mass)/=size(base_mass).or. &
         any(shape(base_velocity)/=[3,size(base_mass)]).or.any(shape(base_velocity)/=shape(velocity_rhs)).or. &
         any(shape(stage_velocity)/=shape(base_velocity))) &
         error stop "block RHS RK extents differ"
    stage_mass=base_mass+h*mass_rhs
    stage_velocity=base_velocity+h*velocity_rhs
  end subroutine
end module parallel_block_rhs_mod
