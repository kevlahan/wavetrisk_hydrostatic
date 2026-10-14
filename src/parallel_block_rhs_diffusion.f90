module parallel_block_rhs_diffusion_mod
  ! Complete level sweeps in compact numerical storage. The halo procedure
  ! binds precompiled producer/consumer plans; it does not interpret arithmetic.
  ! Physical layers are batched because these operators have no vertical reads.
  use kind_mod, only : dp
  use parallel_block_mass_mod, only : Mass_Restriction,execute_mass_restriction
  use parallel_block_rhs_mod, only : Rhs_Primitive_Stencil,Rhs_Pentagon_Stencil
  implicit none
  private
  public :: Rhs_Diffusion_Level,Rhs_Diffusion_Divergence,Rhs_Curl_Stencil
  public :: rhs_laplacian_sweep,rhs_divergence_sweep,rhs_rotational_laplacian
  public :: RHS_DIFFUSION_EDGES,RHS_DIFFUSION_NODES,RHS_ROTATION_EDGES,RHS_ROTATION_TRIANGLES
  public :: RHS_DIFFUSION_DIRECT_EDGES,RHS_ROTATION_DIRECT_TRIANGLES,RHS_DIFFUSION_DIRECT_NODES
  integer, parameter :: RHS_DIFFUSION_DIRECT_EDGES=5,RHS_ROTATION_DIRECT_TRIANGLES=6,RHS_DIFFUSION_DIRECT_NODES=7
  integer, parameter :: RHS_DIFFUSION_EDGES=1,RHS_DIFFUSION_NODES=2,RHS_ROTATION_EDGES=3,RHS_ROTATION_TRIANGLES=4
  integer, parameter :: IJPLUS=5,IPLUSJMINUS=6,IJMINUS=7,IMINUSJPLUS=8
  type :: Rhs_Diffusion_Divergence
     ! Diffusion divergence is UNMASKED. An inactive mass stencil cannot be
     ! reused: its edge addresses may intentionally be uninitialized.
     integer :: target=0,edge(6)=0
     real(dp) :: inverse_area=0.0_dp
  end type
  type :: Rhs_Diffusion_Level
     integer :: level=-1
     type(Rhs_Primitive_Stencil), allocatable :: producer(:)
     type(Mass_Restriction), allocatable :: restriction(:)
     type(Rhs_Diffusion_Divergence), allocatable :: divergence(:)
  end type
  type :: Rhs_Curl_Stencil
     integer :: target=0,vorticity(4)=0
     real(dp) :: dual_length(3)=0.0_dp
     logical :: diagonal=.true.
  end type
  abstract interface
     subroutine Rhs_Diffusion_Halo(phase,level,values,context)
       import dp
       integer, intent(in) :: phase,level
       real(dp), intent(inout) :: values(:,:)
       class(*), intent(inout) :: context
     end subroutine
  end interface
contains
  subroutine rhs_laplacian_sweep(levels,scalar,flux,laplacian,halo,context)
    ! cal_Laplacian_scalars/cal_Laplacian_divu: fine-to-coarse gradient,
    ! flux restriction using completed fine divergence, edge halo, divergence,
    ! then node halo. scalar contains perturbation state or completed divu.
    type(Rhs_Diffusion_Level), intent(in) :: levels(:)
    real(dp), intent(in) :: scalar(:,:)
    real(dp), intent(inout) :: flux(:,:),laplacian(:,:)
    procedure(Rhs_Diffusion_Halo) :: halo
    class(*), intent(inout) :: context
    integer :: l,k
    call check_extents(scalar,flux,laplacian)
    do l=size(levels),1,-1
       do k=1,size(scalar,2)
          call gradient_flux(levels(l)%producer,scalar(:,k),flux(:,k))
       end do
       call halo(RHS_DIFFUSION_DIRECT_EDGES,levels(l)%level,flux,context)
       do k=1,size(scalar,2)
          call execute_mass_restriction(levels(l)%restriction,flux(:,k),laplacian(:,k),record_work=.false.)
       end do
       call halo(RHS_DIFFUSION_EDGES,levels(l)%level,flux,context)
       call divergence(levels(l)%divergence,flux,laplacian)
       call halo(RHS_DIFFUSION_NODES,levels(l)%level,laplacian,context)
    end do
  end subroutine

  subroutine rhs_divergence_sweep(levels,velocity,flux,divu,halo,context)
    ! cal_divu_ml uses the same hierarchy but starts with u*dual_length.
    type(Rhs_Diffusion_Level), intent(in) :: levels(:)
    real(dp), intent(in) :: velocity(:,:)
    real(dp), intent(inout) :: flux(:,:),divu(:,:)
    procedure(Rhs_Diffusion_Halo) :: halo
    class(*), intent(inout) :: context
    integer :: l,k,s,e,n
    if (any(shape(velocity)/=shape(flux))) error stop 'RHS divu velocity extent'
    call check_extents(divu,flux,divu)
    do l=size(levels),1,-1
       do k=1,size(velocity,2)
          do s=1,size(levels(l)%producer)
             associate(p=>levels(l)%producer(s))
             n=3*(p%target-1)
             flux(n+1,k)=velocity(p%velocity(1),k)*p%geometry(11)
             flux(n+2,k)=velocity(p%velocity(4),k)*p%geometry(12)
             flux(n+3,k)=velocity(p%velocity(8),k)*p%geometry(13)
             if (p%southwest) then
                do e=1,3
                   n=3*(p%sw_edge(2,e)-1)+p%sw_edge(1,e)
                   flux(n,k)=velocity(p%velocity(merge(3,merge(6,10,e==2),e==1)),k)*p%geometry(13+e)
                end do
             end if
             end associate
          end do
       end do
       call halo(RHS_DIFFUSION_DIRECT_EDGES,levels(l)%level,flux,context)
       do k=1,size(velocity,2)
          call execute_mass_restriction(levels(l)%restriction,flux(:,k),divu(:,k),record_work=.false.)
       end do
       call halo(RHS_DIFFUSION_EDGES,levels(l)%level,flux,context)
       call divergence(levels(l)%divergence,flux,divu)
       ! No legacy physical node halo is applied. Distributed readers still
       ! require the raw completed fine divergence before coarse restriction.
       call halo(RHS_DIFFUSION_DIRECT_NODES,levels(l)%level,divu,context)
    end do
  end subroutine

  subroutine check_extents(scalar,flux,result)
    real(dp), intent(in) :: scalar(:,:),flux(:,:),result(:,:)
    if (any(shape(scalar)/=shape(result)).or.size(flux,1)/=3*size(scalar,1).or. &
         size(flux,2)/=size(scalar,2)) error stop 'RHS diffusion workspace extent'
  end subroutine

  subroutine gradient_flux(plan,scalar,flux)
    type(Rhs_Primitive_Stencil), intent(in) :: plan(:)
    real(dp), intent(in) :: scalar(:)
    real(dp), intent(inout) :: flux(:)
    real(dp) :: x(4),g(3)
    integer :: s,n,e,edge(3)
    do s=1,size(plan)
       associate(p=>plan(s))
       x=scalar(p%mass([1,3,6,2]))
       g=[x(2)-x(1),x(1)-x(3),x(4)-x(1)]
       n=3*(p%target-1)
       flux(n+1:n+3)=g/p%geometry([1,4,8])*p%geometry(11:13)
       if (p%southwest) then
          x=[scalar(p%mass(1)),scalar(p%mass(5)),scalar(p%sw_mass),scalar(p%mass(4))]
          g=[x(2)-x(1),x(1)-x(3),x(4)-x(1)]
          edge=3*(p%sw_edge(2,:)-1)+p%sw_edge(1,:)
          ! Keep the accepted SW sign order (gradient negation before length).
          g=-g/p%geometry([3,6,10])*p%geometry(14:16)
          do e=1,3
             flux(edge(e))=g(e)
          end do
       end if
       end associate
    end do
  end subroutine

  subroutine divergence(plan,flux,result)
    type(Rhs_Diffusion_Divergence), intent(in) :: plan(:)
    real(dp), intent(in) :: flux(:,:)
    real(dp), intent(inout) :: result(:,:)
    integer :: k,s
    do k=1,size(flux,2)
       do s=1,size(plan)
          associate(p=>plan(s))
          result(p%target,k)=(flux(p%edge(1),k)-flux(p%edge(2),k)+flux(p%edge(3),k)-flux(p%edge(4),k) &
               +flux(p%edge(5),k)-flux(p%edge(6),k))*p%inverse_area
          end associate
       end do
    end do
  end subroutine

  subroutine rhs_rotational_laplacian(level,curl_plan,producer,pentagon,vorticity,rotation,halo,context)
    ! cal_Laplacian_vector_rot: curl(vort), edge halo, ordered circulation
    ! including SW writes and post_vort pentagons, then final-writer delivery.
    integer, intent(in) :: level
    type(Rhs_Curl_Stencil), intent(in) :: curl_plan(:)
    type(Rhs_Primitive_Stencil), intent(in) :: producer(:)
    type(Rhs_Pentagon_Stencil), intent(in) :: pentagon(:)
    real(dp), intent(inout) :: vorticity(:,:),rotation(:,:)
    procedure(Rhs_Diffusion_Halo) :: halo
    class(*), intent(inout) :: context
    real(dp) :: x(4),u(10),sw(2)
    integer :: k,s,n
    if (size(vorticity,2)/=size(rotation,2).or.3*size(vorticity,1)/=2*size(rotation,1)) &
         error stop 'RHS rotation workspace extent'
    do k=1,size(vorticity,2)
       do s=1,size(curl_plan)
          associate(p=>curl_plan(s))
          x=vorticity(p%vorticity,k);n=3*(p%target-1)
          rotation(n+1,k)=-(x(1)-x(4))/p%dual_length(1)
          rotation(n+2,k)=0.0_dp
          if (p%diagonal) rotation(n+2,k)=-(x(1)-x(2))/p%dual_length(2)
          rotation(n+3,k)=-(x(3)-x(2))/p%dual_length(3)
          end associate
       end do
    end do
    call halo(RHS_ROTATION_EDGES,level,rotation,context)
    do k=1,size(vorticity,2)
       do s=1,size(producer)
          associate(p=>producer(s))
          u=rotation(p%velocity,k)*p%geometry(1:10)
          vorticity(p%vorticity(1),k)=(u(1)+u(9)+u(4))/p%geometry(34)
          vorticity(p%vorticity(2),k)=-(u(4)+u(8)+u(2))/p%geometry(35)
          vorticity(p%vorticity(3),k)=(u(3)+u(8)+u(7))/p%geometry(36)
          vorticity(p%vorticity(4),k)=-(u(1)+u(5)+u(10))/p%geometry(37)
          if (p%southwest) then
             sw=rotation(p%sw_velocity,k)*p%sw_geometry(1:2)
             vorticity(p%sw_vorticity(1),k)=(sw(1)+u(10)+u(6))/p%sw_geometry(11)
             vorticity(p%sw_vorticity(2),k)=-(u(3)+u(6)+sw(2))/p%sw_geometry(12)
          end if
          end associate
       end do
    end do
    call halo(RHS_ROTATION_DIRECT_TRIANGLES,level,vorticity,context)
    do k=1,size(vorticity,2)
       do s=1,size(pentagon)
          call pentagon_rotation(pentagon(s),rotation(:,k),vorticity(:,k))
       end do
    end do
    call halo(RHS_ROTATION_TRIANGLES,level,vorticity,context)
  end subroutine

  subroutine pentagon_rotation(p,velocity,vorticity)
    type(Rhs_Pentagon_Stencil), intent(in) :: p
    real(dp), intent(in) :: velocity(:)
    real(dp), intent(inout) :: vorticity(:)
    real(dp) :: u(3,7),x
    integer :: n
    do n=1,7
       u(:,n)=velocity(p%velocity(:,n))*p%length(:,n)
    end do
    select case(p%corner)
    case(IJMINUS)
       x=(u(3,4)-u(1,5)+u(1,7))/p%triangle_area(1,7)
       vorticity(p%vorticity(1,7))=x;vorticity(p%vorticity(2,7))=x
    case(IPLUSJMINUS)
       x=(-u(1,1)+u(1,7)+u(2,7))/p%triangle_area(1,7)
       vorticity(p%vorticity(1,7))=x;vorticity(p%vorticity(2,4))=x
    case(IMINUSJPLUS)
       x=(u(3,1)-u(2,7)-u(3,7))/p%triangle_area(2,7)
       vorticity(p%vorticity(2,7))=x;vorticity(p%vorticity(1,5))=x
    case(IJPLUS)
       x=(u(1,1)-u(1,2)-u(3,1))/p%triangle_area(1,1)
       vorticity(p%vorticity(1,1))=x;vorticity(p%vorticity(2,1))=x
    case default
       error stop 'RHS rotational pentagon invalid'
    end select
  end subroutine
end module
