  ! Validation-only prototype: direct fused arithmetic on resident block state.
  ! Reference packet traffic is an oracle, never a proposed production exchange.
  integer, parameter :: P276_WIDTH=68
  integer, parameter :: p276_dx(7)=[0,0,1,0,-1,1,-1],p276_dy(7)=[0,1,0,-1,0,1,-1]
  integer, parameter :: p276_vnode(10)=[1,2,5,1,4,7,5,1,3,4]
  integer, parameter :: p276_vedge(10)=[1,1,1,2,2,2,2,3,3,3]
  integer, parameter :: p276_dual(6)=[1,4,8,3,6,10]
  integer, parameter :: p276_inode(16)=[1,2,3,4,5,6,p276_vnode],p276_iedge(16)=[0,0,0,0,0,0,p276_vedge]
  type :: P276_Domain_Map
     integer, allocatable :: b(:),p(:)
  end type
  type :: Primitive276_Plan
     integer, allocatable :: storage(:,:),node(:,:)
     real(dp), allocatable :: output(:,:)
  end type
  type(P276_Domain_Map), allocatable, save :: p276_map(:)
  type(Primitive276_Plan), allocatable, save :: p276_plan(:)
  integer(int64), save :: p276_generation=-1_int64
  integer, save :: p276_transaction=0,p276_unit=-1,p276_used=0
  logical, save :: p276_active=.false.
  real(dp), allocatable, save :: p276_send(:,:),p276_recv(:,:)
  integer, allocatable, save :: p276_dest(:)
  real(dp), save :: p276_prepare=0,p276_capture_time=0,p276_exchange=0,p276_address=0, &
       p276_gather=0,p276_compute=0,p276_check=0
  integer(int64), save :: p276_inputs(3)=0,p276_diff(3)=0,p276_outputs=0,p276_output_diff=0, &
       p276_ref_diff=0,p276_nodes=0,p276_skipped=0
