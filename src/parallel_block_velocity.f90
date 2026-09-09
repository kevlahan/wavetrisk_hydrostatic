module parallel_block_velocity_mod
  ! Native velocity operators and persistent numerical workspace. No Domain
  ! or MPI dependency: adapters supply primitives at their restriction phase.
  use, intrinsic :: iso_fortran_env, only : int64
  use kind_mod, only : dp
  implicit none
  private
  public :: velocity_dx, velocity_dy, velocity_weights, velocity_qperp
  public :: velocity_source, velocity_restrict_source, velocity_gradient
  public :: velocity_pressure_layer
  public :: Velocity_Direct_Stencil, Velocity_Restriction_Stencil
  public :: Velocity_Gradient_Stencil, Velocity_Level_Program
  public :: VELOCITY_DIRECT, VELOCITY_RESTRICT
  public :: validate_velocity_program, execute_velocity_sources, execute_velocity_gradients
  public :: Native_Velocity_Workspace, native_velocity, native_velocity_work
  public :: native_velocity_transaction, native_velocity_oracle
  public :: copy_native_velocity_tendency, native_velocity_rk

  integer, parameter :: VELOCITY_DIRECT = 1, VELOCITY_RESTRICT = 2

  ! Addresses are one-based indices in the caller's native workspace, not
  ! Domain identifiers. A route compiler may remap them to any execution owner.
  type :: Velocity_Direct_Stencil
     integer :: node(9) = 0
     real(dp) :: weights(5,2,3) = 0.0_dp
     real(dp) :: length(3) = 0.0_dp
     logical :: active = .false.
  end type

  type :: Velocity_Restriction_Stencil
     integer :: target = 0, child = 0, neighbor(3) = 0
     logical :: edge(3) = .false.
  end type

  type :: Velocity_Gradient_Stencil
     integer :: node(4) = 0
     real(dp) :: length(3) = 0.0_dp
     logical :: active = .false.
  end type

  type :: Velocity_Level_Program
     ! Ordered action tape: direct writes and restricted writes must not be
     ! regrouped into separate passes. Mixed masks may require a direct write
     ! immediately before a partial restriction at the same target.
     integer, allocatable :: action(:), operand(:)
     type(Velocity_Direct_Stencil), allocatable :: direct(:)
     type(Velocity_Restriction_Stencil), allocatable :: restriction(:)
  end type

  type :: Native_Velocity_Workspace
     integer(int64) :: generation = -1_int64
     type(Velocity_Level_Program), allocatable :: level(:)
     type(Velocity_Gradient_Stencil), allocatable :: gradient(:)
     integer, allocatable :: physics_address(:,:)
     real(dp), allocatable :: source(:,:), physics(:,:), tendency(:,:), rk(:,:), rho(:), rho_theta(:)
     logical, allocatable :: physics_ready(:), ready(:)
     logical, allocatable :: active_mask(:), restriction_mask(:)
  end type
  type(Native_Velocity_Workspace), allocatable, target, save :: native_velocity(:)
  integer(int64), save :: native_velocity_work(9) = 0_int64
  logical, save :: native_velocity_transaction = .false., native_velocity_oracle = .false.

  ! Center, E, NE, N, W, SW, S, SE, NW. Edges are RT, DG, UP.
  integer, parameter :: velocity_dx(9) = [0,1,1,0,-1,-1,0,1,-1]
  integer, parameter :: velocity_dy(9) = [0,0,1,1,0,-1,-1,-1,1]

contains

  subroutine copy_native_velocity_tendency(d,k,first,values)
    integer, intent(in) :: d,k,first
    real(dp), intent(out) :: values(:)
    call assert_native_velocity_ready(d,k,first,size(values))
    values = native_velocity(d)%tendency(first:first+size(values)-1,k)
    native_velocity_work(6) = native_velocity_work(6)+int(size(values),int64)
  end subroutine

  subroutine native_velocity_rk(d,k,first,solution,h,values)
    integer, intent(in) :: d,k,first
    real(dp), intent(in) :: solution(:),h
    real(dp), intent(out) :: values(:)
    call assert_native_velocity_ready(d,k,first,size(values))
    if (size(solution) /= size(values)) error stop "native velocity RK extent differs"
    values = solution+h*native_velocity(d)%tendency(first:first+size(values)-1,k)
    native_velocity_work(7) = native_velocity_work(7)+int(size(values),int64)
  end subroutine

  subroutine assert_native_velocity_ready(d,k,first,n)
    integer, intent(in) :: d,k,first,n
    if (.not. allocated(native_velocity)) error stop "native velocity workspace is absent"
    if (d < 1 .or. d > size(native_velocity)) error stop "native velocity workspace owner is invalid"
    if (.not. allocated(native_velocity(d)%ready)) error stop "native velocity readiness is absent"
    if (k < 1 .or. k > size(native_velocity(d)%ready)) error stop "native velocity layer is invalid"
    if (.not. native_velocity(d)%ready(k)) error stop "native velocity tendency is not ready"
    if (first < 1 .or. n < 0 .or. first+n-1 > size(native_velocity(d)%tendency,1)) &
         error stop "native velocity tendency extent is invalid"
  end subroutine

  subroutine validate_velocity_program(program,nnode)
    type(Velocity_Level_Program), intent(in) :: program
    integer, intent(in) :: nnode
    integer :: a,s,e
    if (nnode < 0) error stop "velocity program: negative workspace extent"
    if (.not. allocated(program%action) .or. .not. allocated(program%operand) .or. &
         .not. allocated(program%direct) .or. .not. allocated(program%restriction)) &
         error stop "velocity program: incomplete plan"
    if (size(program%action) /= size(program%operand)) error stop "velocity program: action extent differs"
    do a = 1,size(program%action)
       s = program%operand(a)
       select case(program%action(a))
       case(VELOCITY_DIRECT)
          if (s < 1 .or. s > size(program%direct)) error stop "velocity program: invalid direct action"
          associate(stencil=>program%direct(s))
          if (stencil%node(1) < 1 .or. stencil%node(1) > nnode) error stop "velocity program: invalid direct target"
          if (stencil%active) then
             if (any(stencil%node < 1) .or. any(stencil%node > nnode)) &
                  error stop "velocity program: incomplete primitive neighborhood"
             if (any(stencil%length <= 0.0_dp)) error stop "velocity program: nonpositive direct metric"
          end if
          end associate
       case(VELOCITY_RESTRICT)
          if (s < 1 .or. s > size(program%restriction)) error stop "velocity program: invalid restriction action"
          associate(stencil=>program%restriction(s))
          if (stencil%target < 1 .or. stencil%target > nnode) error stop "velocity program: invalid restriction target"
          do e = 1,3
             if (.not. stencil%edge(e)) cycle
             if (stencil%child < 1 .or. stencil%child > nnode) error stop "velocity program: invalid child"
             if (stencil%neighbor(e) < 1 .or. stencil%neighbor(e) > nnode) &
                  error stop "velocity program: invalid child neighbor"
          end do
          end associate
       case default
          error stop "velocity program: unknown action"
       end select
    end do
  end subroutine validate_velocity_program

  subroutine execute_velocity_sources(program,flux,pv,physics,source,qperp_values)
    type(Velocity_Level_Program), intent(in) :: program
    real(dp), intent(in) :: flux(:,:),pv(:,:),physics(:,:)
    real(dp), intent(inout) :: source(:,:)
    real(dp), optional, intent(inout) :: qperp_values(:,:)
    real(dp) :: qperp(3)
    integer :: a,s,target,e,nnode
    nnode = size(source,2)
    if (size(source,1) /= 3 .or. any(shape(flux) /= [3,nnode]) .or. &
         any(shape(pv) /= [3,nnode]) .or. any(shape(physics) /= [3,nnode])) &
         error stop "velocity source execution: workspace extents differ"
    if (present(qperp_values)) then
       if (any(shape(qperp_values) /= [3,nnode])) error stop "velocity source execution: oracle extent differs"
    end if
    do a = 1,size(program%action)
       s = program%operand(a)
       select case(program%action(a))
       case(VELOCITY_DIRECT)
          associate(stencil=>program%direct(s))
          target = stencil%node(1)
          qperp = 0.0_dp
          if (stencil%active) then
             qperp = velocity_qperp(flux(:,stencil%node),pv(:,stencil%node),stencil%weights)
             source(:,target) = velocity_source(qperp,physics(:,target),stencil%length,.true.)
          else
             source(:,target) = 0.0_dp
          end if
          if (present(qperp_values)) qperp_values(:,target) = qperp
          end associate
       case(VELOCITY_RESTRICT)
          associate(stencil=>program%restriction(s))
          do e = 1,3
             if (stencil%edge(e)) &
                  source(e,stencil%target) = source(e,stencil%child)+source(e,stencil%neighbor(e))
          end do
          end associate
       case default
          error stop "velocity source execution: unknown action"
       end select
    end do
  end subroutine execute_velocity_sources

  subroutine execute_velocity_gradients(plan,bernoulli,exner,mass,temperature,source,tendency)
    type(Velocity_Gradient_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: bernoulli(:),exner(:),mass(:),temperature(:),source(:,:)
    real(dp), intent(inout) :: tendency(:,:)
    integer :: s,target,nnode
    nnode = size(source,2)
    if (size(source,1) /= 3 .or. any(shape(tendency) /= [3,nnode]) .or. &
         size(bernoulli) /= nnode .or. size(exner) /= nnode .or. &
         size(mass) /= nnode .or. size(temperature) /= nnode) &
         error stop "velocity gradient execution: workspace extents differ"
    do s = 1,size(plan)
       target = plan(s)%node(1)
       if (plan(s)%active) then
          tendency(:,target) = velocity_gradient(source(:,target),plan(s)%length, &
               bernoulli(plan(s)%node),exner(plan(s)%node),mass(plan(s)%node),temperature(plan(s)%node),.true.)
       else
          tendency(:,target) = 0.0_dp
       end if
    end do
  end subroutine execute_velocity_gradients

  pure function velocity_weights(parts,inverse_area) result(weights)
    ! Compile once per topology generation. The four columns are C,E,NE,N.
    real(dp), intent(in) :: parts(6,4), inverse_area(4)
    real(dp) :: weights(5,2,3), cumulative(5)
    integer :: edge, side, node, offset, term
    do edge = 1,3
       do side = 1,2
          node = 1
          offset = edge-1
          if (side == 2) then
             node = edge+1
             offset = edge+2
          end if
          cumulative(1) = parts(1+offset,node)
          do term = 2,5
             cumulative(term) = cumulative(term-1) + parts(modulo(term+offset-1,6)+1,node)
          end do
          cumulative = 0.5_dp-cumulative*inverse_area(node)
          weights(:,side,edge) = [cumulative(1),-cumulative(2),cumulative(3),-cumulative(4),cumulative(5)]
       end do
    end do
  end function velocity_weights

  pure function velocity_qperp(flux,pv,weights) result(value)
    real(dp), intent(in) :: flux(3,9), pv(3,9), weights(5,2,3)
    real(dp) :: value(3)
    ! Keep the ten-term order in ops_mod:Qperp, including multiplication by
    ! the average before the weight. Do not use a reduction/reassociated sum.
    value(1) = &
         flux(2,1)*(0.5_dp*(pv(2,1)+pv(1,1)))*weights(1,1,1) + &
         flux(3,1)*(0.5_dp*(pv(3,1)+pv(1,1)))*weights(2,1,1) + &
         flux(1,5)*(0.5_dp*(pv(1,5)+pv(1,1)))*weights(3,1,1) + &
         flux(2,6)*(0.5_dp*(pv(2,6)+pv(1,1)))*weights(4,1,1) + &
         flux(3,7)*(0.5_dp*(pv(3,7)+pv(1,1)))*weights(5,1,1) + &
         flux(2,7)*(0.5_dp*(pv(2,7)+pv(1,1)))*weights(1,2,1) + &
         flux(3,8)*(0.5_dp*(pv(3,8)+pv(1,1)))*weights(2,2,1) + &
         flux(1,2)*(0.5_dp*(pv(1,2)+pv(1,1)))*weights(3,2,1) + &
         flux(2,2)*(0.5_dp*(pv(2,2)+pv(1,1)))*weights(4,2,1) + &
         flux(3,2)*(0.5_dp*(pv(3,2)+pv(1,1)))*weights(5,2,1)
    value(2) = &
         flux(3,1)*(0.5_dp*(pv(3,1)+pv(2,1)))*weights(1,1,2) + &
         flux(1,5)*(0.5_dp*(pv(1,5)+pv(2,1)))*weights(2,1,2) + &
         flux(2,6)*(0.5_dp*(pv(2,6)+pv(2,1)))*weights(3,1,2) + &
         flux(3,7)*(0.5_dp*(pv(3,7)+pv(2,1)))*weights(4,1,2) + &
         flux(1,1)*(0.5_dp*(pv(1,1)+pv(2,1)))*weights(5,1,2) + &
         flux(3,2)*(0.5_dp*(pv(3,2)+pv(2,1)))*weights(1,2,2) + &
         flux(1,3)*(0.5_dp*(pv(1,3)+pv(2,1)))*weights(2,2,2) + &
         flux(2,3)*(0.5_dp*(pv(2,3)+pv(2,1)))*weights(3,2,2) + &
         flux(3,3)*(0.5_dp*(pv(3,3)+pv(2,1)))*weights(4,2,2) + &
         flux(1,4)*(0.5_dp*(pv(1,4)+pv(2,1)))*weights(5,2,2)
    value(3) = &
         flux(1,5)*(0.5_dp*(pv(1,5)+pv(3,1)))*weights(1,1,3) + &
         flux(2,6)*(0.5_dp*(pv(2,6)+pv(3,1)))*weights(2,1,3) + &
         flux(3,7)*(0.5_dp*(pv(3,7)+pv(3,1)))*weights(3,1,3) + &
         flux(1,1)*(0.5_dp*(pv(1,1)+pv(3,1)))*weights(4,1,3) + &
         flux(2,1)*(0.5_dp*(pv(2,1)+pv(3,1)))*weights(5,1,3) + &
         flux(1,4)*(0.5_dp*(pv(1,4)+pv(3,1)))*weights(1,2,3) + &
         flux(2,4)*(0.5_dp*(pv(2,4)+pv(3,1)))*weights(2,2,3) + &
         flux(3,4)*(0.5_dp*(pv(3,4)+pv(3,1)))*weights(3,2,3) + &
         flux(1,9)*(0.5_dp*(pv(1,9)+pv(3,1)))*weights(4,2,3) + &
         flux(2,5)*(0.5_dp*(pv(2,5)+pv(3,1)))*weights(5,2,3)
  end function velocity_qperp

  pure function velocity_source(qperp,physics,length,active) result(value)
    real(dp), intent(in) :: qperp(3), physics(3), length(3)
    logical, intent(in) :: active
    real(dp) :: value(3)
    value = 0.0_dp
    if (active) value = -qperp+physics*length
  end function velocity_source

  pure function velocity_restrict_source(direct,child,neighbor,restrict_edge) result(value)
    ! Direct evaluation must occur for an absent child or any child mask
    ! below ADJZONE. The caller schedules those cases; this kernel selects
    ! each edge independently, without reading an unused child component.
    real(dp), intent(in) :: direct(3), child(3), neighbor(3)
    logical, intent(in) :: restrict_edge(3)
    real(dp) :: value(3)
    integer :: edge
    do edge = 1,3
       if (restrict_edge(edge)) then
          value(edge) = child(edge)+neighbor(edge)
       else
          value(edge) = direct(edge)
       end if
    end do
  end function velocity_restrict_source

  pure function velocity_gradient(source,length,bernoulli,exner,mass,temperature,active) result(value)
    ! Scalar columns are C,E,NE,N. B/Exner must be the FINAL restricted fields,
    ! not an Exner field reconstructed independently at the consumer.
    real(dp), intent(in) :: source(3), length(3), bernoulli(4), exner(4), mass(4), temperature(4)
    logical, intent(in) :: active
    real(dp) :: value(3), theta(4), theta_edge(3), grad_b(3), grad_e(3)
    value = 0.0_dp
    if (.not. active) return
    theta = temperature/mass
    theta_edge = 0.5_dp*(theta(1)+theta(2:4))
    grad_b = [bernoulli(2)-bernoulli(1),bernoulli(1)-bernoulli(3),bernoulli(4)-bernoulli(1)]/length
    grad_e = [exner(2)-exner(1),exner(1)-exner(3),exner(4)-exner(1)]/length
    value = source/length-grad_b-theta_edge*grad_e
  end function velocity_gradient

  pure subroutine velocity_pressure_layer(mass,temperature,gravity,cp,p0,kappa, &
       lower_pressure,lower_geopotential,pressure,upper_pressure,exner,upper_geopotential)
    ! Compressible physical-layer producer. Surface pressure and surface
    ! geopotential are inputs; the caller advances layers in physical order.
    ! Preserve integrate_pressure_up's midpoint, not lower - 0.5*g*mass.
    real(dp), intent(in) :: mass, temperature, gravity, cp, p0, kappa, lower_pressure, lower_geopotential
    real(dp), intent(out) :: pressure, upper_pressure, exner, upper_geopotential
    upper_pressure = lower_pressure-gravity*mass
    pressure = 0.5_dp*(lower_pressure+upper_pressure)
    exner = cp*(pressure/p0)**kappa
    upper_geopotential = lower_geopotential+gravity*kappa*temperature*exner/pressure
  end subroutine velocity_pressure_layer

end module parallel_block_velocity_mod
