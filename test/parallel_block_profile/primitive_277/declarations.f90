  ! Stage277: cached numerical operands only; no ownership or transport change.
  integer, parameter :: P277_DX(7)=[0,0,1,0,-1,1,-1],P277_DY(7)=[0,1,0,-1,0,1,-1]
  integer, parameter :: P277_VNODE(10)=[1,2,5,1,4,7,5,1,3,4],P277_VEDGE(10)=[1,1,1,2,2,2,2,3,3,3]
  integer, parameter :: P277_DUAL(6)=[1,4,8,3,6,10]
  type :: Primitive277_Patch
     integer, allocatable :: velocity(:,:),mass(:,:)
     real(dp), allocatable :: geometry(:,:)
  end type
  type :: Primitive277_Domain
     type(Primitive277_Patch), allocatable :: patch(:)
  end type
  type(Primitive277_Domain), allocatable, save :: primitive277(:)
  integer(int64), save :: p277_generation=-1_int64,p277_rebuilds=0,p277_checks=0
  integer(int64), allocatable, save :: p277_layer_nodes(:)
  integer, save :: p277_transaction=0,p277_unit=-1
  logical, save :: p277_initialized=.false.,p277_check=.false.
