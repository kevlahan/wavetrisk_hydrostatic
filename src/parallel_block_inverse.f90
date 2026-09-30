module parallel_block_inverse_mod
  use parallel_block_profile_mod
  ! Native inverse boundary dependency workspace. Topology is compiled once
  ! per generation; only required nodes are staged, never Float_Field copies.
  ! Regular scalar/inner-vector kernels remain on final-owner blocks. The
  ! sparse outer-edge producer retains the established geometric route owner
  ! and write order, including aliases which have no compact interior owner.
  use iso_fortran_env, only : int64, error_unit
  use ieee_arithmetic, only : ieee_is_finite
  use mpi_f08, only : MPI_Alltoall, MPI_Alltoallv, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUCCESS, MPI_Wtime, &
       MPI_Request, MPI_Irecv, MPI_Isend, MPI_Waitall, MPI_STATUSES_IGNORE, MPI_Allreduce, MPI_MIN, MPI_MAX
  use kind_mod, only : dp
  use shared_mod, only : EDGE, AT_NODE, AT_EDGE, N_VARIABLE, N_BDRY, N_CHDRN, &
       scalars, S_VELO, S_TEMP, zlevels, zmin, n_domain, level_start, level_end, end_pt, opp_no, hex_sides, hex_s_offs, &
       ADJZONE, Z_NULL, DG, FROZEN, LORT, UPLT, TRIAG, RT, UP, NORTHEAST, NORTHWEST, IMINUSJPLUS, IPLUSJMINUS, IJMINUS
  use domain_mod, only : grid, sol, sol_mean, Float_Field, idx, ed_idx, get_offs_Domain, chd_offs, idx__fast
  use patch_mod, only : PATCH_SIZE, LAST
  use arch_mod, only : comm, rank, n_process, glo_id, loc_id, owner, block_catalog, abort_run
  use parallel_block_mod, only : Block_Data, STORE_PATCH, STORE_BDRY, BLOCK_PAYLOAD_SOL, &
       BLOCK_PAYLOAD_WAV_COEFF, apply_local_block_field_consumer, &
       compile_local_block_inverse_routes, transfer_local_block_inverse_routes
  use parallel_block_primitive_walk_mod, only: walk_primitive_patch
  use domain_ops_mod, only: apply_onescale_to_patch,apply_interscale_d,apply_to_penta_d,apply_onescale_d,apply_interscale_d2
  use domain_mod, only: Domain
  implicit none
  private
  public :: prepare_native_inverse, native_inverse_gather, native_inverse_boundary, &
       native_inverse_scatter, native_inverse_outer, compare_native_inverse, publish_native_inverse
  public :: inverse_operand_key, inverse_scalar_operand, inverse_vector_operand, inverse_operand_generation
  public :: complete_inverse_operand_epoch, require_provisional_inverse_operands
  public :: require_tendency_operands,prepare_initial_tendency_operands
  public :: capture_native_inverse_consumer_aliases
  public :: physics_consumer_alias,finish_physics_consumer_epoch,install_physics_consumer_alias
  public :: export_native_soil_aliases,native_soil_transform,publish_native_soil
  real(dp), public :: native_inverse_seconds(3)=0.0_dp
  integer(int64), public :: native_inverse_calls(3)=0_int64, native_inverse_messages(3)=0_int64, &
       native_inverse_bytes(3)=0_int64
#ifdef BLOCK_RK
  logical, parameter, public :: block_rk_enabled=.true.
#else
  logical, parameter, public :: block_rk_enabled=.false.
#endif
#if defined(BLOCK_RK) && !defined(BLOCK_COLUMNS_ORACLE) && !defined(BLOCK_PRIMITIVES_ORACLE) && !defined(BLOCK_DIFFUSION_ORACLE) && !defined(BLOCK_MASS_OWNER_ORACLE) && !defined(BLOCK_VELOCITY_OWNER_ORACLE)
  logical,parameter,public::block_rk_native=.true.
#else
  logical,parameter,public::block_rk_native=.false.
#endif
  public :: prepare_retained_rk_aliases,finish_retained_rk_boundary,rk_scaffold_operand
  public :: reset_retained_rk_scaffold
  integer,allocatable::rk_low_key(:,:),rk_low_slot(:)
  real(dp),allocatable::rk_low_value(:,:,:)
  integer::original_inverse_nodes=0
  public :: retained_rk_cycle_active
  public :: compare_retained_rk_aliases
  public :: rk_alias_nodes,set_rk_mass_alias_layer,rk_mass_alias_operand
  real(dp), allocatable :: rk_alias_base(:,:,:),rk_mass_rhs(:,:)
  integer,allocatable::rk_alias_index(:)
  logical,allocatable::rk_mass_ready(:)
  integer :: rk_stage=0,rk_stages=0
  logical :: rk_boundary_ready=.false.,rk_base_ready=.false.
  logical :: provisional_operands_ready=.false.,initial_operands_ready=.false.
  logical :: timing=.false.
  logical :: sparse_outer_install=.false., stage192_ready=.false.

  type :: Node_Map
     integer, allocatable :: slot(:), key(:,:)
  end type
  type :: Transfer_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),key(:,:),slot(:)
     integer, allocatable :: source_info(:,:),target_info(:,:)
     integer, allocatable :: address(:,:)
#ifdef STAGE202_INVERSE
     integer, allocatable :: self_slot(:)
#endif
     logical :: edgewise=.false.
  end type
  type :: Alias_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),source(:),dest(:)
     integer, allocatable :: install_order(:)
  end type
  type :: Local_Aliases
     integer, allocatable :: source(:),dest(:)
  end type
  type :: Outer_Operation
     integer :: level=0,kind=0,source(13)=0,target(3)=0
     real(dp) :: weight(9)=0.0_dp
  end type
  type :: Boundary_Manifest
     integer, allocatable :: item(:,:)
  end type
  type(Node_Map), allocatable :: nodes(:)
  type(Transfer_Plan) :: interior_plan,boundary_plan
  type(Transfer_Plan), allocatable :: level_interior(:,:),level_boundary(:,:)
  type(Transfer_Plan), allocatable :: outer_install(:)
  type(Alias_Plan), allocatable :: aliases(:,:)
  type(Alias_Plan) :: range_alias(2,2)
  type(Local_Aliases) :: local_alias(2)
  type(Outer_Operation), allocatable :: operations(:)
  type(MPI_Request), allocatable :: requests(:)
  real(dp), allocatable, asynchronous :: transfer_send(:),transfer_recv(:)
  real(dp), allocatable, asynchronous :: alias_send(:),alias_recv(:)
  integer, allocatable :: node_domain(:),node_id(:),operation_first(:),operation_last(:)
  integer, allocatable :: scaffold_slot(:)
  integer, allocatable :: consumer_alias_slot(:),consumer_alias_index(:)
  logical :: physics_consumer_ready=.false.
  real(dp), allocatable :: consumer_alias_value(:,:,:)
  real(dp), allocatable :: physics_alias_mean(:,:,:),physics_alias_soil(:,:,:)
  logical :: physics_reference_ready=.false.
  integer, allocatable :: node_flags(:)
  real(dp), allocatable :: scaffold_value(:,:,:)
  integer(int64), allocatable :: coverage(:,:)
  real(dp), allocatable,target :: value(:,:,:,:)
  integer :: soil_domain=0,soil_layer=0,soil_variable=0
  integer :: inverse_field_first=1,inverse_levels=0
  logical :: physics_results_ready=.false.
  integer(int64) :: plan_generation=-1_int64
  integer :: nnode=0,noperation=0,nscalar=0
  ! The selector is compile-time; the original path remains the default.
#ifdef STAGE203_INVERSE
  logical, parameter, public :: stage203_inverse=.true.
#else
  logical, parameter, public :: stage203_inverse=.false.
#endif
  public :: begin_inverse203, end_inverse203, scalar_inverse203, inner_inverse203, install_inverse203
#ifdef STAGE203_LEDGER
  logical, parameter, public :: ledger203=.true.
#else
  logical, parameter, public :: ledger203=.false.
#endif
  logical :: resident203=.false., oracle203=.false.
  type :: Scalar203
     integer :: d=0,level=0,phase=0,target=0,base=0,area=0,mask=0,n=0
     integer :: source(11)=0,geometry(11)=0
  end type
  type :: Inner203
     integer :: d=0,level=0,node(32)=0,edge(32)=0,source(32)=0,target(6)=0
  end type
  type(Scalar203), allocatable :: scalar203_ops(:)
  type(Inner203), allocatable :: inner203_ops(:)
  integer, allocatable :: restore203_index(:),restore203_slot(:)
  real(dp), allocatable :: restore203(:,:,:),saved203(:,:,:)
  integer(int64), public :: work203(6)=0_int64
  real(dp), parameter :: base_weight(9) = &
       [16.0_dp,-1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,-1.0_dp,1.0_dp,1.0_dp]/16.0_dp

contains
  subroutine reset_retained_rk_scaffold
    if(allocated(rk_low_key)) deallocate(rk_low_key,rk_low_value)
    if(allocated(rk_low_slot)) deallocate(rk_low_slot)
    rk_base_ready=.false.;rk_boundary_ready=.false.;rk_stage=0;rk_stages=0
  end subroutine
  subroutine prepare_rk_low_scaffold
    integer,allocatable::keys(:,:)
    real(dp),allocatable::values(:,:,:)
    integer::q,p,d,id,n,t,v,key(2)
    if(allocated(rk_low_slot)) deallocate(rk_low_slot)
    allocate(rk_low_slot(nnode));rk_low_slot=0
    n=0
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q);d=node_domain(p);id=node_id(p)
       if(id<grid(d)%patch%elts(3)%elts_start) n=n+1
    end do
    allocate(keys(2,n),values(nscalar+EDGE,zlevels,n));n=0
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q);d=node_domain(p);id=node_id(p)
       if(id>=grid(d)%patch%elts(3)%elts_start) cycle
       n=n+1;key=[glo_id(rank+1,d),id];keys(:,n)=key;rk_low_slot(p)=n
       do v=1,nscalar
          values(v,:,n)=real(N_VARIABLE-(scalars(1)+v-1),dp)
       end do
       values(nscalar+1:,:,n)=real(N_VARIABLE-S_VELO,dp)
       if(.not.allocated(rk_low_key)) cycle
       do t=1,size(rk_low_key,2)
          if(all(rk_low_key(:,t)==key)) then
             values(:,:,n)=rk_low_value(:,:,t)
             exit
          end if
       end do
    end do
    call move_alloc(keys,rk_low_key);call move_alloc(values,rk_low_value)
  end subroutine

  logical function retained_rk_cycle_active() result(active)
    active=block_rk_enabled.and.rk_base_ready.and.rk_stage>0
  end function

  subroutine rk_alias_nodes(address)
    integer,allocatable,intent(out)::address(:,:)
    integer::q,p
    allocate(address(2,size(consumer_alias_slot)))
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q);address(:,q)=[node_domain(p),node_id(p)]
    end do
  end subroutine
  subroutine set_rk_mass_alias_layer(k,rhs)
    integer,intent(in)::k
    real(dp),intent(in)::rhs(:)
    if(.not.rk_base_ready.or.k<1.or.k>zlevels) call die('RK mass alias epoch')
    if(rk_mass_ready(k).or.size(rhs)/=size(rk_mass_rhs,2)) call die('RK mass alias duplicate or extent')
    rk_mass_rhs(k,:)=rhs;rk_mass_ready(k)=.true.
  end subroutine
  real(dp) function rk_mass_alias_operand(d,id,k) result(x)
    integer,intent(in)::d,id,k
    integer::p,q
    if(.not.rk_base_ready.or..not.rk_mass_ready(k)) call die('RK mass alias layer unavailable')
    p=nodes(d)%slot(id+1)
    if(p==0) call die('RK mass alias not registered')
    q=rk_alias_index(p)
    if(q==0) call die('RK mass alias requested canonical interior')
    x=rk_mass_rhs(k,q)
  end function

  subroutine prepare_retained_rk_aliases(h,stage,stages,rhs,vector_closure)
    real(dp),intent(in)::h
    integer,intent(in)::stage,stages
    interface
       subroutine rhs(d,id,k,tendency)
         import dp
         integer,intent(in)::d,id,k
         real(dp),intent(out)::tendency(:)
       end subroutine
       function vector_closure(key,k,h) result(x)
         import dp
         integer,intent(in)::key(4),k
         real(dp),intent(in)::h
         real(dp)::x(3)
       end function
    end interface
    integer::q,p,d,id,k,ierr,n
    real(dp)::tendency(nscalar+EDGE)
    if(.not.block_rk_enabled.or..not.rk_base_ready) call die('RK base is unavailable')
    if(rk_boundary_ready.or.stage/=rk_stage+1) call die('RK closure stage sequence')
    if(stages/=3.and.stages/=4) call die('RK closure method')
    if(stage<1.or.stage>stages) call die('RK closure stage range')
    if(stage>1.and.stages/=rk_stages) call die('RK closure method changed')
    rk_stage=stage;rk_stages=stages
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q);d=node_domain(p);id=node_id(p)
       do k=1,zlevels
          call rhs(d,id,k,tendency)
          if(.not.all(ieee_is_finite(tendency))) call die('RK alias missing tendency')
          value(:,k,p,BLOCK_PAYLOAD_SOL)=rk_alias_base(:,k,q)+h*tendency
          if(stage<stages.and.rk_low_slot(p)>0) then
             value(:,k,p,BLOCK_PAYLOAD_SOL)=rk_low_value(:,k,rk_low_slot(p))
             ! Native temperature publication preserves global-sol scaffold,
             ! unlike the untouched provisional mass/velocity storage.
             if(nodes(d)%key(2,id+1)==-1) &
                  value(S_TEMP-scalars(1)+1,k,p,BLOCK_PAYLOAD_SOL)=rk_alias_base(S_TEMP-scalars(1)+1,k,q)
          end if
       end do
    end do
    if(stage==stages) then
       do q=1,size(scaffold_slot)
          scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),BLOCK_PAYLOAD_SOL)
       end do
    end if
    call native_inverse_gather(1,BLOCK_PAYLOAD_SOL,.false.)
    ! The canonical block candidate includes the accepted add/subtract
    ! recomposition. Physical aliases require the original RK rounding order.
    ! Read the existing owner checkpoint and exact retained gradient directly.
    n=EDGE*zlevels
    do p=1,size(interior_plan%key,2)
       do k=1,zlevels
          q=(p-1)*n+(k-1)*EDGE
          transfer_recv(q+1:q+EDGE)=vector_closure(interior_plan%key(:,p),k,h)
       end do
    end do
    call MPI_Alltoallv(transfer_recv,n*interior_plan%rc,n*interior_plan%rd,MPI_DOUBLE_PRECISION, &
         transfer_send,n*interior_plan%sc,n*interior_plan%sd,MPI_DOUBLE_PRECISION,comm,ierr)
    call check(ierr,'RK exact vector closure')
    do p=1,size(interior_plan%slot)
       do k=1,zlevels
          q=(p-1)*n+(k-1)*EDGE
          value(nscalar+1:nscalar+EDGE,k,interior_plan%slot(p),BLOCK_PAYLOAD_SOL)=transfer_send(q+1:q+EDGE)
       end do
    end do
    if(stage==stages) then
       do d=1,size(grid)
          do k=1,zlevels
             call apply_interscale_d(restrict_rk_velocity,grid(d),level_start-1,k,0,0)
          end do
       end do
       do q=1,size(scaffold_slot)
          scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),BLOCK_PAYLOAD_SOL)
       end do
    end if
    rk_boundary_ready=.true.
    rk_mass_ready=.false.
    write(*,'(a,i0,a,i0,a,2(1x,i0))') 'Retained RK closure rank=',rank,' stage=',stage, &
         ' aliases/owned operands=',size(consumer_alias_slot),size(interior_plan%key,2)
  contains
    subroutine restrict_rk_velocity(dom,ip,jp,ic,jc,k,op,dp_,oc,dc)
      type(Domain),intent(inout)::dom
      integer,intent(in)::ip,jp,ic,jc,k,op(N_BDRY+1),dp_(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
      integer::target,c,e,ne,n,ts,cs,es,nes,ns
      target=idx(ip,jp,op,dp_);c=idx(ic,jc,oc,dc)
      e=idx(ic+1,jc,oc,dc);ne=idx(ic+1,jc+1,oc,dc);n=idx(ic,jc+1,oc,dc)
      ts=nodes(d)%slot(target+1);cs=nodes(d)%slot(c+1);es=nodes(d)%slot(e+1)
      nes=nodes(d)%slot(ne+1);ns=nodes(d)%slot(n+1)
      if(min(ts,cs,es,nes,ns)<1) call die('RK coarse restriction operand missing')
      if(dom%mask_e%elts(EDGE*c+RT+1)>0) value(nscalar+RT+1,k,ts,1)= &
           0.5_dp*(value(nscalar+RT+1,k,cs,1)+value(nscalar+RT+1,k,es,1))
      if(dom%mask_e%elts(EDGE*c+DG+1)>0) value(nscalar+DG+1,k,ts,1)= &
           0.5_dp*(value(nscalar+DG+1,k,nes,1)+value(nscalar+DG+1,k,cs,1))
      if(dom%mask_e%elts(EDGE*c+UP+1)>0) value(nscalar+UP+1,k,ts,1)= &
           0.5_dp*(value(nscalar+UP+1,k,cs,1)+value(nscalar+UP+1,k,ns,1))
    end subroutine
  end subroutine

  subroutine compare_retained_rk_aliases(reference)
    type(Float_Field),intent(in)::reference(1:N_VARIABLE,1:zlevels)
    integer::q,p,d,id,k,v,var,address
    real(dp)::expected,actual
    if(.not.rk_boundary_ready) call die('RK comparison outside boundary epoch')
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q)
       if(p>original_inverse_nodes) cycle
       d=node_domain(p);id=node_id(p)
       do k=1,zlevels
          do v=1,nscalar+EDGE
             if(v<=nscalar) then
                var=scalars(1)+v-1;address=id+1
             else
                var=S_VELO;address=EDGE*id+v-nscalar
             end if
             expected=reference(var,k)%data(d)%elts(address)
             actual=value(v,k,p,BLOCK_PAYLOAD_SOL)
             if(transfer_bits(actual)==transfer_bits(expected)) cycle
             write(error_unit,*) 'RK alias exact mismatch rank/stage/d/id/k/v:',rank,rk_stage,d,id,k,v
             write(error_unit,*) 'actual/reference:',actual,expected
             call die('independent RK alias mismatch')
          end do
       end do
    end do
    write(*,'(a,i0,a,i0,a,i0)') 'Retained RK exact aliases rank=',rank,' stage=',rk_stage,' count=', &
         size(consumer_alias_slot)*zlevels*(nscalar+EDGE)
  end subroutine

  subroutine finish_retained_rk_boundary
    integer::q
    if(.not.rk_boundary_ready) call die('RK boundary before closure')
    call native_inverse_boundary(1,BLOCK_PAYLOAD_SOL,level_start-1,level_end)
    call native_inverse_boundary(2,BLOCK_PAYLOAD_SOL,level_start-1,level_end)
    call native_inverse_scatter(1,BLOCK_PAYLOAD_SOL,.false.)
    call native_inverse_scatter(2,BLOCK_PAYLOAD_SOL,.false.)
    if(rk_stage==rk_stages) then
       ! Final closure replaces global sol, including signed boundary aliases.
       ! Subsequent scaffold gathers must use that completed boundary epoch.
       do q=1,size(scaffold_slot)
          scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),BLOCK_PAYLOAD_SOL)
       end do
    end if
  end subroutine

  real(dp) function rk_scaffold_operand(d,id,component,k) result(x)
    integer,intent(in)::d,id,component,k
    integer::p,q
    if(.not.rk_base_ready) call die('RK scaffold without base')
    p=nodes(d)%slot(id+1)
    if(p==0) call die('RK scaffold was not registered')
    if(nodes(d)%key(1,id+1)/=0) call die('RK scaffold requested canonical interior')
    do q=1,size(scaffold_slot)
       if(scaffold_slot(q)==p) then
          x=scaffold_value(component,k,q)
          return
       end if
    end do
    call die('RK scaffold absent')
  end function

  subroutine complete_inverse_operand_epoch(provisional)
    logical, intent(in) :: provisional
    integer::p
    if(block_rk_enabled.and.provisional) then
       do p=1,original_inverse_nodes
          if(rk_low_slot(p)>0) rk_low_value(:,:,rk_low_slot(p))=value(:,:,p,BLOCK_PAYLOAD_SOL)
       end do
    end if
    initial_operands_ready=.false.
    provisional_operands_ready=provisional
    if(.not.provisional) rk_base_ready=.false.
  end subroutine

  subroutine require_provisional_inverse_operands
    if(.not.provisional_operands_ready) call die('provisional consumer operand epoch is not ready')
  end subroutine

  subroutine require_tendency_operands
    if(inverse_field_first/=1.or.inverse_levels/=zlevels) call die('non-atmospheric tendency epoch')
    if(.not.provisional_operands_ready.and..not.initial_operands_ready) call die('tendency operand epoch is not ready')
  end subroutine

  subroutine prepare_initial_tendency_operands(generation,physical)
    ! Explicit input boundary for stage one: the upstream restart/adaptation/
    ! physics producer still exports physical aliases through its compatibility
    ! field. Only registered aliases/scaffold enter this epoch. Canonical
    ! interiors always come from final-owner Block_Data. This does NOT select
    ! the distinct saved final-inverse consumer epoch or perform another inverse.
    integer(int64),intent(in)::generation
    type(Float_Field),intent(in)::physical(1:N_VARIABLE,1:zlevels)
    call prepare_native_inverse(generation,physical,physical,level_start,.false.,consumer_only=.true.)
  end subroutine

  function inverse_operand_generation() result(generation)
    integer(int64) :: generation
    generation=plan_generation
  end function

  subroutine inverse_operand_key(d,id,key)
    integer, intent(in) :: d,id
    integer, intent(out) :: key(4)
    if (.not.allocated(nodes)) call die('operand plan is absent')
    if (d<1.or.d>size(nodes)) call die('operand Domain address')
    if (id<0.or.id>=size(nodes(d)%slot)) call die('operand node address')
    key=nodes(d)%key(:,id+1)
  end subroutine

  function inverse_vector_operand(d,id,k) result(operand)
    integer, intent(in) :: d,id,k
    real(dp) :: operand(EDGE)
    integer :: key(4),slot
    call require_tendency_operands
    call inverse_operand_key(d,id,key)
    if(initial_operands_ready.and.key(1)>0) call die('initial epoch attempted an interior compatibility read')
    slot=nodes(d)%slot(id+1)
    if(slot==0) call die('required vector operand was not registered')
    if(k<1.or.k>zlevels) call die('vector operand layer')
    operand=value(nscalar+1:nscalar+EDGE,k,slot,BLOCK_PAYLOAD_SOL)
    if(.not.all(ieee_is_finite(operand))) call die('nonfinite resident vector operand')
  end function

  real(dp) function inverse_scalar_operand(d,id,v,k) result(operand)
    ! Named provisional or initial-physical epoch; never the saved final-
    ! inverse consumer alias version.
    integer, intent(in) :: d,id,v,k
    integer :: key(4),slot
    call require_tendency_operands
    call inverse_operand_key(d,id,key)
    if(initial_operands_ready.and.key(1)>0) call die('initial epoch attempted an interior compatibility read')
    slot=nodes(d)%slot(id+1)
    if(slot==0) call die('required scalar operand was not registered')
    if(v<scalars(1).or.v>scalars(2).or.k<1.or.k>zlevels) call die('scalar operand component')
    operand=value(v-scalars(1)+1,k,slot,BLOCK_PAYLOAD_SOL)
    if(.not.ieee_is_finite(operand)) call die('nonfinite resident scalar operand')
  end function

#ifdef BLOCK_COLUMNS
  subroutine register_column_operand(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: slot
    slot=node_slot(dom%id+1,idx(i,j,offs,dims))
  end subroutine
#endif

  subroutine configure_stage192
    character(16) :: setting
    integer :: status,length,selected,low,high,ierr
    if (stage192_ready) return
    call get_environment_variable('WAVETRISK_STAGE192_SPARSE_OUTER',setting,length,status)
    if (status/=1) then
       if (status/=0 .or. length/=1) call die('STAGE192_SPARSE_OUTER must be 0 or 1')
       select case(setting(1:1))
       case('0')
          sparse_outer_install=.false.
       case('1')
          sparse_outer_install=.true.
       case default
          call die('STAGE192_SPARSE_OUTER must be 0 or 1')
       end select
    end if
    selected=merge(1,0,sparse_outer_install)
    call MPI_Allreduce(selected,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call check(ierr,'Stage 192 selector minimum')
    call MPI_Allreduce(selected,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call check(ierr,'Stage 192 selector maximum')
    if (low/=high) call die('Stage 192 selector differs between ranks')
    if (rank==0) write(*,'(a,i0)') 'Stage 192 requested sparse outer install = ',selected
    stage192_ready=.true.
  end subroutine

  subroutine check(ierr,where)
    integer, intent(in) :: ierr
    character(*), intent(in) :: where
    if (ierr /= MPI_SUCCESS) call die(where)
  end subroutine

  subroutine die(message)
    character(*), intent(in) :: message
    write(error_unit,'(a,i0,2a)') 'Rank ',rank,': native inverse: ',message
    call abort_run
    error stop "native inverse: MPI abort returned"
  end subroutine

  subroutine displacements(count,offset)
    integer, intent(in) :: count(:)
    integer, intent(out) :: offset(:)
    integer :: r
    offset=0
    do r=2,size(count)
       offset(r)=offset(r-1)+count(r-1)
    end do
  end subroutine

  integer function node_slot(d,id) result(slot)
    integer, intent(in) :: d,id
    if (d < 1 .or. d > size(nodes)) call die('invalid geometry owner')
    if (id < 0 .or. id >= size(nodes(d)%slot)) call die('invalid geometry node')
    slot=nodes(d)%slot(id+1)
    if (slot > 0) return
    nnode=nnode+1
    slot=nnode
    nodes(d)%slot(id+1)=slot
    node_domain(slot)=d
    node_id(slot)=id
  end function

  integer function edge_slot(d,id) result(slot)
    integer, intent(in) :: d,id
    integer :: q
    if (id < 0) call die('negative source edge')
    q=node_slot(d,id/EDGE)
    slot=EDGE*(q-1)+mod(id,EDGE)+1
  end function

  integer function route_slot(d,id,kind) result(slot)
    integer, intent(in) :: d,id,kind
    if (kind == 1) then
       slot=node_slot(d,abs(id))
    else
       slot=edge_slot(d,abs(id))
       if (id < 0) slot=-slot
    end if
  end function

  subroutine prepare_native_inverse(generation,scaling,wavelet,first_level,profile,consumer_only, &
       field_first,field_count,entry_alias_address,entry_alias_value)
    integer(int64), intent(in) :: generation
    integer, intent(in) :: first_level
    logical, intent(in) :: profile
    logical,intent(in),optional :: consumer_only
    logical :: initial_consumer
    type(Float_Field), intent(in),optional :: scaling(:,:),wavelet(:,:)
    integer,intent(in),optional::field_first,field_count,entry_alias_address(:,:)
    real(dp),intent(in),optional::entry_alias_value(:,:,:)
    integer::selected_first,selected_count,old_first,old_count
    logical::native_entry
    logical,allocatable::entry_seen(:)
    integer :: d,id,q,k,v,total,b,p,r,pos,n,next_patch
    integer, allocatable :: sc(:),rc(:),sd(:),rd(:),sb(:,:),rb(:,:),request(:,:),slots(:)
    type(Boundary_Manifest) :: manifest
    type(Transfer_Plan) :: outer_patch,outer_bdry

    selected_first=1;selected_count=zlevels
    if(present(field_first)) selected_first=field_first
    if(present(field_count)) selected_count=field_count
    if(present(field_first).neqv.present(field_count)) call die('incomplete inverse vertical range')
    if(selected_count<1.or.selected_first<zmin.or.selected_first+selected_count-1>zlevels) &
         call die('invalid inverse vertical range')
    native_entry=present(entry_alias_value)
    if(native_entry.neqv.present(entry_alias_address)) call die('incomplete native inverse entry')
    if(.not.native_entry) then
       if(.not.present(scaling).or..not.present(wavelet)) call die('inverse compatibility inputs absent')
       if(selected_first/=1.or.selected_count/=zlevels) call die('non-atmospheric inverse requires native entry')
       if(size(scaling,1)/=N_VARIABLE.or.size(scaling,2)/=selected_count.or. &
            size(wavelet,1)/=N_VARIABLE.or.size(wavelet,2)/=selected_count) call die('inverse input shape')
    end if
    old_first=inverse_field_first;old_count=inverse_levels
    inverse_field_first=selected_first;inverse_levels=selected_count
    physics_results_ready=.false.
    provisional_operands_ready=.false.;initial_operands_ready=.false.
    physics_consumer_ready=.false.
    initial_consumer=.false.
    if(present(consumer_only)) initial_consumer=consumer_only
    resident203=.false.
    call configure_stage192
    call detail_enter(DP_INVERSE_PLAN)
    timing=profile
    native_inverse_seconds=0.0_dp
    native_inverse_calls=0_int64
    native_inverse_messages=0_int64
    native_inverse_bytes=0_int64
    if(initial_consumer.and.(native_entry.or.selected_first/=1.or.selected_count/=zlevels)) &
         call die('initial tendency requires atmospheric entry')
    if (plan_generation /= generation.or.old_first/=selected_first.or.old_count/=selected_count) then
       call detail_add(DC_INVERSE_BUILD,1_int64)
       if (allocated(nodes)) deallocate(nodes,node_domain,node_id,value,aliases,operations,coverage, &
            scaffold_slot,scaffold_value,operation_first,operation_last,node_flags,level_interior,level_boundary, &
            outer_install,transfer_send,transfer_recv,alias_send,alias_recv,consumer_alias_slot,consumer_alias_value, &
            consumer_alias_index,physics_alias_mean,physics_alias_soil)
       do q=1,2
          if (allocated(local_alias(q)%source)) deallocate(local_alias(q)%source,local_alias(q)%dest)
       end do
       allocate(nodes(size(grid)))
       if (.not. allocated(requests)) allocate(requests(2*n_process))
       total=0
       nscalar=scalars(2)-scalars(1)+1
       do d=1,size(grid)
          n=grid(d)%node%length
          allocate(nodes(d)%slot(n),nodes(d)%key(4,n))
          nodes(d)%slot=0
          nodes(d)%key=0
          do p=1,grid(d)%patch%length
             if (grid(d)%patch%elts(p)%deleted) cycle
             id=grid(d)%patch%elts(p)%elts_start
             nodes(d)%key(2,id+1:id+PATCH_SIZE**2)=-1
          end do
          total=total+n
       end do
       allocate(node_domain(total),node_id(total))
       nnode=0
       do b=1,size(block_catalog)
          if (owner(block_catalog(b)%root_domain+1) /= rank) cycle
          d=loc_id(block_catalog(b)%root_domain+1)+1
          next_patch=0
          call map_patch(block_catalog(b)%root_patch)
       end do
       if (stage203_inverse) call compile203(total)
       allocate(aliases(2,level_start-1:level_end))
       call compile_aliases
       do q=1,2
          call combine_aliases(range_alias(q,1),q,level_start-1)
          call combine_aliases(range_alias(q,2),q,level_start)
       end do
       allocate(alias_send(max(1,inverse_levels*max(nscalar*sum(range_alias(1,1)%sc),sum(range_alias(2,1)%sc)))))
       allocate(alias_recv(max(1,inverse_levels*max(nscalar*sum(range_alias(1,1)%rc),sum(range_alias(2,1)%rc)))))
       ! The upper bound is conservative; records are reused for all fields.
       allocate(operations(max(1,EDGE*total)))
       allocate(coverage(4,level_start-1:level_end))
       allocate(operation_first(level_start-1:level_end),operation_last(level_start-1:level_end))
       operation_first=1
       operation_last=0
       coverage=0_int64
       noperation=0
       call compile_outer

       ! Final owners supply the compact boundary destinations. Only integer
       ! geometry metadata uses original Domain indices; values never do.
       allocate(manifest%item(6,0))
       call apply_local_block_field_consumer(collect_boundaries,manifest)
       allocate(sc(n_process),rc(n_process),sd(n_process),rd(n_process))
       sc=0
       do p=1,size(manifest%item,2)
          r=owner(manifest%item(5,p)+1)+1
          sc(r)=sc(r)+1
       end do
       call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,q)
       call check(q,'boundary manifest counts')
       call displacements(sc,sd)
       call displacements(rc,rd)
       allocate(sb(6,max(1,sum(sc))),rb(6,max(1,sum(rc))))
       sc=0
       do p=1,size(manifest%item,2)
          r=owner(manifest%item(5,p)+1)+1
          pos=sd(r)+sc(r)+1
          sb(:,pos)=manifest%item(:,p)
          sc(r)=sc(r)+1
       end do
       call MPI_Alltoallv(sb,6*sc,6*sd,MPI_INTEGER,rb,6*rc,6*rd,MPI_INTEGER,comm,q)
       call check(q,'boundary manifest')
       n=sum(rc)
       allocate(request(4,n),slots(n))
       do p=1,n
          d=loc_id(rb(5,p)+1)+1
          slots(p)=node_slot(d,rb(6,p))
          request(:,p)=rb(1:4,p)
       end do
       ! Flags are carried in the existing request-key exchange. Both ends
       ! then derive identical per-level subplans without extra collectives.
#ifdef BLOCK_COLUMNS
       ! Register the exact pressure consumer closure, including repeated aliases.
       ! No all-geometry solution snapshot is created.
       do d=1,size(grid)
          do p=3,grid(d)%patch%length
             call apply_onescale_to_patch(register_column_operand,grid(d),p-1,0,0,1)
          end do
       end do
#endif
#ifdef BLOCK_PRIMITIVES
       do d=1,size(grid)
          do p=3,grid(d)%patch%length
             if(grid(d)%patch%elts(p)%deleted) cycle
             call walk_primitive_patch(grid(d),p-1,register_primitive_operand)
          end do
       end do
#endif
       original_inverse_nodes=nnode
       if(block_rk_enabled) then
          ! Fixed coarse patches remain explicit scaffold consumers of the
          ! forward transform. Register this bounded subset, never interiors.
          do d=1,size(nodes)
             do id=0,size(nodes(d)%slot)-1
                if(nodes(d)%key(2,id+1)==-1) r=node_slot(d,id)
             end do
          end do
       end if
       allocate(node_flags(nnode))
       node_flags=0
       do q=1,2
          do p=1,size(local_alias(q)%dest)
             r=alias_node(local_alias(q)%dest(p),q)
             node_flags(r)=ibset(node_flags(r),q-1)
             node_flags(r)=ibset(node_flags(r),q+2)
             if (q==2) node_flags(r)=ibset(node_flags(r),6+mod(abs(local_alias(q)%dest(p))-1,EDGE))
          end do
          do b=level_start-1,level_end
             do p=1,size(aliases(q,b)%dest)
                r=alias_node(aliases(q,b)%dest(p),q)
                node_flags(r)=ibset(node_flags(r),q-1)
                if (q==2) node_flags(r)=ibset(node_flags(r),6+mod(abs(aliases(q,b)%dest(p))-1,EDGE))
             end do
          end do
       end do
       do p=1,nnode
          if (nodes(node_domain(p))%key(2,node_id(p)+1)==-1) node_flags(p)=ibset(node_flags(p),2)
       end do
       ! Plus-side or pentagon targets outside the nominal fine level must
       ! participate too; do not assume geometry level alone is sufficient.
       do p=1,noperation
          do q=2,3
             r=alias_node(operations(p)%target(q),2)
             ! Only target(2:3) are written by native_inverse_outer. Keep a
             ! conservative union across levels, shared with final owners in
             ! the existing key exchange; target(1) is a read-only parent.
             node_flags(r)=ibset(node_flags(r),5)
             node_flags(r)=ibset(node_flags(r),6+mod(abs(operations(p)%target(q))-1,EDGE))
             if (grid(node_domain(r))%level%elts(node_id(r)+1)/=operations(p)%level+1) then
                node_flags(r)=ibset(node_flags(r),1)
                node_flags(r)=ibset(node_flags(r),4)
             end if
          end do
       end do
       call build_transfer(boundary_plan,request,slots)
       deallocate(request,slots)
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(1,node_id(p)+1) > 0) n=n+1
       end do
       allocate(request(4,n),slots(n))
       n=0
       do p=1,nnode
          d=node_domain(p)
          id=node_id(p)
          if (nodes(d)%key(1,id+1) == 0) cycle
          n=n+1
          request(:,n)=nodes(d)%key(:,id+1)
          slots(n)=p
       end do
       call build_transfer(interior_plan,request,slots)
       ! Every transfer completes before the next one starts. Subplans keep
       ! only integer routing metadata, sharing this one maximum-size pair.
       allocate(transfer_send(max(1,max(nscalar,EDGE)*inverse_levels*(sum(interior_plan%sc)+sum(boundary_plan%sc)))))
       allocate(transfer_recv(max(1,max(nscalar,EDGE)*inverse_levels*(sum(interior_plan%rc)+sum(boundary_plan%rc)))))
       allocate(level_interior(2,level_start-2:level_end),level_boundary(2,level_start-2:level_end))
       do q=1,2
          do b=level_start-2,level_end
             call select_transfer(level_interior(q,b),interior_plan,q,b,.false.,.false.)
             call select_transfer(level_boundary(q,b),boundary_plan,q,b,.true.,.false.)
          end do
       end do
       allocate(outer_install(level_start-1:level_end))
       do b=level_start-1,level_end
          call select_transfer(outer_patch,interior_plan,2,b,.true.,.true.,sparse_outer_install)
          call select_transfer(outer_bdry,boundary_plan,2,b,.true.,.true.,sparse_outer_install)
          if (sparse_outer_install) then
             call combine_outer_edges(outer_install(b),outer_patch,outer_bdry)
          else
             call combine_transfers(outer_install(b),outer_patch,outer_bdry)
          end if
       end do
#ifdef STAGE202_INVERSE
       call prepare_self_routes(interior_plan)
       call prepare_self_routes(boundary_plan)
       do q=1,2
          do b=level_start-2,level_end
             call prepare_self_routes(level_interior(q,b))
             call prepare_self_routes(level_boundary(q,b))
          end do
       end do
       do b=level_start-1,level_end
          call prepare_self_routes(outer_install(b))
       end do
#endif
       allocate(value(nscalar+EDGE,inverse_levels,nnode,2))
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(1,node_id(p)+1)==0) n=n+1
       end do
       allocate(consumer_alias_slot(n),consumer_alias_value(nscalar+EDGE,inverse_levels,n),consumer_alias_index(nnode))
#ifdef PHYSICS_NATIVE_STATE
       allocate(physics_alias_mean(nscalar,inverse_levels,n))
#ifdef SOIL_NATIVE
       allocate(physics_alias_soil(nscalar+EDGE,1-zmin,n))
#else
       allocate(physics_alias_soil(S_TEMP-scalars(1)+1:S_TEMP-scalars(1)+1,1-zmin,n))
#endif
#else
       allocate(physics_alias_mean(0,0,0),physics_alias_soil(0,0,0))
#endif
       physics_reference_ready=.false.
       consumer_alias_index=0
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(1,node_id(p)+1)/=0) cycle
          n=n+1
          consumer_alias_slot(n)=p
          consumer_alias_index(p)=n
       end do
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(2,node_id(p)+1)==-1) n=n+1
       end do
       allocate(scaffold_slot(n),scaffold_value(nscalar+EDGE,inverse_levels,n))
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(2,node_id(p)+1)/=-1) cycle
          n=n+1
          scaffold_slot(n)=p
       end do
       if (stage203_inverse) then
          if (allocated(restore203_index)) deallocate(restore203_index,restore203_slot,restore203,saved203)
          allocate(restore203_index(nnode));restore203_index=0
          n=0
          do p=1,nnode
             if (nodes(node_domain(p))%key(1,node_id(p)+1)<=0) cycle
             if (.not.btest(node_flags(p),0).and..not.btest(node_flags(p),1)) cycle
             n=n+1;restore203_index(p)=n
          end do
          allocate(restore203_slot(n),restore203(nscalar+EDGE,inverse_levels,n),saved203(nscalar+EDGE,inverse_levels,n))
          do p=1,nnode
             if (restore203_index(p)>0) restore203_slot(restore203_index(p))=p
          end do
       end if
       plan_generation=generation
    else
       call detail_add(DC_INVERSE_REUSE,1_int64)
    end if
    call detail_leave(DP_INVERSE_PLAN)
    if(native_entry) then
       if(size(entry_alias_address,1)/=2.or.size(entry_alias_value,1)/=nscalar+EDGE.or. &
            size(entry_alias_value,2)/=inverse_levels.or. &
            size(entry_alias_value,3)/=size(entry_alias_address,2)) call die('native inverse entry shape')
       value=0.0_dp
       allocate(entry_seen(nnode));entry_seen=.false.
       do q=1,size(entry_alias_address,2)
          d=entry_alias_address(1,q);id=entry_alias_address(2,q)
          if(d<1.or.d>size(nodes)) call die('native inverse entry Domain')
          if(id<0.or.id>=size(nodes(d)%slot)) call die('native inverse entry node')
          p=nodes(d)%slot(id+1)
          if(p<1.or.nodes(d)%key(1,id+1)>0) call die('native inverse entry is not an alias')
          if(entry_seen(p)) call die('duplicate native inverse entry')
          value(:,:,p,BLOCK_PAYLOAD_SOL)=entry_alias_value(:,:,q);entry_seen(p)=.true.
       end do
       do q=1,size(consumer_alias_slot)
          p=consumer_alias_slot(q)
          if(.not.entry_seen(p)) call die('native inverse alias entry missing')
       end do
       do q=1,size(scaffold_slot)
          scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),BLOCK_PAYLOAD_SOL)
       end do
       do q=1,2
          call native_inverse_gather(q,BLOCK_PAYLOAD_SOL,.false.)
       end do
       if(stage203_inverse) then
          do q=1,size(restore203_slot)
             restore203(:,:,q)=value(:,:,restore203_slot(q),BLOCK_PAYLOAD_SOL)
          end do
       end if
       do q=1,size(consumer_alias_slot)
          consumer_alias_value(:,:,q)=value(:,:,consumer_alias_slot(q),BLOCK_PAYLOAD_SOL)
       end do
       return
    end if
    if(initial_consumer) then
       if(block_rk_enabled) then
          if(rk_boundary_ready) call die('RK boundary not consumed before next step')
          if(allocated(rk_alias_base)) deallocate(rk_alias_base,rk_mass_rhs,rk_mass_ready,rk_alias_index)
          allocate(rk_mass_rhs(inverse_levels,size(consumer_alias_slot)),rk_mass_ready(inverse_levels),rk_alias_index(nnode))
          rk_mass_ready=.false.;rk_alias_index=0
          do q=1,size(consumer_alias_slot)
             rk_alias_index(consumer_alias_slot(q))=q
          end do
          allocate(rk_alias_base(nscalar+EDGE,inverse_levels,size(consumer_alias_slot)))
          rk_stage=0;rk_stages=0;rk_base_ready=.true.
          call prepare_rk_low_scaffold
       end if
       do q=1,size(consumer_alias_slot)
          p=consumer_alias_slot(q);d=node_domain(p);id=node_id(p)
          if(nodes(d)%key(1,id+1)/=0) call die('initial alias includes canonical interior')
#ifdef PHYSICS_NATIVE_STATE
          do k=zmin,0
#ifdef SOIL_NATIVE
             do v=1,nscalar
                physics_alias_soil(v,1-k,q)=sol(scalars(1)+v-1,k)%data(d)%elts(id+1)
             end do
             physics_alias_soil(nscalar+1:nscalar+EDGE,1-k,q)= &
                  sol(S_VELO,k)%data(d)%elts(EDGE*id+1:EDGE*id+EDGE)
#else
             physics_alias_soil(S_TEMP-scalars(1)+1,1-k,q)=sol(S_TEMP,k)%data(d)%elts(id+1)
#endif
          end do
          do k=1,inverse_levels
             do v=1,nscalar
                physics_alias_mean(v,k,q)=sol_mean(scalars(1)+v-1,k)%data(d)%elts(id+1)
             end do
          end do
#endif
          do k=1,inverse_levels
             do v=1,nscalar
                value(v,k,p,BLOCK_PAYLOAD_SOL)=scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)
             end do
             value(nscalar+1:nscalar+EDGE,k,p,BLOCK_PAYLOAD_SOL)=scaling(S_VELO,k)%data(d)%elts(EDGE*id+1:EDGE*id+EDGE)
          end do
       end do
       if(block_rk_enabled) then
          do q=1,size(consumer_alias_slot)
             rk_alias_base(:,:,q)=value(:,:,consumer_alias_slot(q),BLOCK_PAYLOAD_SOL)
          end do
          do q=1,size(scaffold_slot)
             scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),BLOCK_PAYLOAD_SOL)
          end do
       end if
       initial_operands_ready=.true.
#ifdef PHYSICS_NATIVE_STATE
       physics_reference_ready=.true.
#endif
       write(*,'(a,i0,a,i0)') 'Initial tendency alias epoch rank=',rank,' registered aliases=',size(consumer_alias_slot)
       return
    end if
    call detail_enter(DP_INVERSE_SEED)

    ! Seed only compatibility/scaffold and pre-existing boundary state once
    ! per transform. Every catalogued interior dependency is replaced below
    ! from the authoritative final owner before any arithmetic or alias copy.
    do p=1,nnode
       d=node_domain(p)
       id=node_id(p)
       if (nodes(d)%key(1,id+1)>0) cycle
       do k=1,inverse_levels
          do v=1,nscalar
             if(.not.block_rk_enabled) value(v,k,p,1)=scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)
             value(v,k,p,2)=wavelet(scalars(1)+v-1,k)%data(d)%elts(id+1)
          end do
          do v=1,EDGE
             if(.not.block_rk_enabled) value(nscalar+v,k,p,1)=scaling(S_VELO,k)%data(d)%elts(EDGE*id+v)
             value(nscalar+v,k,p,2)=wavelet(S_VELO,k)%data(d)%elts(EDGE*id+v)
          end do
       end do
    end do
    ! Legacy writeback preserves uncatalogued interior patches from global
    ! sol, NOT from the provisional scaling argument. Keep this fixed input
    ! separate from the pre-boundary seed and apply it at each gather point.
    if(.not.block_rk_enabled) then
    do q=1,size(scaffold_slot)
       p=scaffold_slot(q)
       d=node_domain(p)
       id=node_id(p)
       do k=1,inverse_levels
          do v=1,nscalar
             scaffold_value(v,k,q)=sol(scalars(1)+v-1,k)%data(d)%elts(id+1)
          end do
          do v=1,EDGE
             scaffold_value(nscalar+v,k,q)=sol(S_VELO,k)%data(d)%elts(EDGE*id+v)
          end do
       end do
    end do
    end if
    if(block_rk_enabled) then
       if(.not.rk_boundary_ready) call die('inverse requires retained RK boundary')
       rk_boundary_ready=.false.
    end if
    do q=1,2
       do b=1,2
          call native_inverse_gather(q,b,.false.)
       end do
    end do
    if (stage203_inverse) then
       do q=1,size(restore203_slot)
          restore203(:,:,q)=value(:,:,restore203_slot(q),1)
       end do
    end if
    do q=1,2
       call native_inverse_boundary(q,2,level_start,level_end)
       call native_inverse_boundary(q,1,first_level,level_end)
       call native_inverse_scatter(q,2,.false.)
       call native_inverse_scatter(q,1,.false.)
    end do
    ! Preserve the legacy consumer-visible alias epoch separately from the
    ! fresher aliases needed by native kernels. Only boundary/scaffold values
    ! are retained, not a second solution or another inverse computation.
    do q=1,size(consumer_alias_slot)
       consumer_alias_value(:,:,q)=value(:,:,consumer_alias_slot(q),BLOCK_PAYLOAD_SOL)
    end do

    call detail_leave(DP_INVERSE_SEED)

  contains
#ifdef BLOCK_PRIMITIVES
    subroutine register_primitive_operand(ids,southwest)
      integer, intent(in) :: ids(0:7)
      logical, intent(in) :: southwest
      integer :: q,slot
      do q=0,7
         slot=node_slot(d,ids(q))
      end do
    end subroutine
#endif
    recursive subroutine map_patch(patch)
      integer, intent(in) :: patch
      integer :: c,child,i,j,id_local
      if (grid(d)%patch%elts(patch+1)%deleted) return
      next_patch=next_patch+1
      do j=0,PATCH_SIZE-1
         do i=0,PATCH_SIZE-1
            id_local=grid(d)%patch%elts(patch+1)%elts_start+PATCH_SIZE*j+i
            nodes(d)%key(:,id_local+1)=[b,STORE_PATCH,next_patch,PATCH_SIZE*j+i]
         end do
      end do
      do c=1,N_CHDRN
         child=grid(d)%patch%elts(patch+1)%children(c)
         if (child > 0) call map_patch(child)
      end do
    end subroutine
  end subroutine

  subroutine collect_boundaries(catalog,block,context)
    integer, intent(in) :: catalog
    type(Block_Data), intent(in) :: block
    class(*), intent(inout) :: context
    integer :: p,q,n,pos
    integer, allocatable :: grown(:,:)
    select type(context)
    type is(Boundary_Manifest)
       n=size(context%item,2)
       pos=sum(block%bdry_storage%n_node)
       allocate(grown(6,n+pos))
       grown(:,1:n)=context%item
       do p=1,size(block%bdry_storage)
          do q=0,block%bdry_storage(p)%n_node-1
             n=n+1
             grown(:,n)=[catalog,STORE_BDRY,p,q,block%root_domain,block%bdry_storage(p)%elts_start+q]
          end do
       end do
       call move_alloc(grown,context%item)
    class default
       call die('boundary manifest context')
    end select
  end subroutine

  subroutine build_transfer(plan,request,slots)
    type(Transfer_Plan), intent(out) :: plan
    integer, intent(in) :: request(:,:),slots(:)
    integer :: p,r,pos,ierr
    integer, allocatable :: send_key(:,:),recv_key(:,:)
    allocate(plan%sc(n_process),plan%sd(n_process),plan%rc(n_process),plan%rd(n_process))
    plan%sc=0
    do p=1,size(slots)
       r=block_catalog(request(1,p))%owner+1
       plan%sc(r)=plan%sc(r)+1
    end do
    call MPI_Alltoall(plan%sc,1,MPI_INTEGER,plan%rc,1,MPI_INTEGER,comm,ierr)
    call check(ierr,'native dependency counts')
    call displacements(plan%sc,plan%sd)
    call displacements(plan%rc,plan%rd)
    allocate(send_key(6,max(1,sum(plan%sc))),recv_key(6,max(1,sum(plan%rc))))
    allocate(plan%key(4,sum(plan%rc)),plan%source_info(2,sum(plan%sc)),plan%target_info(2,sum(plan%rc)))
    allocate(plan%slot(sum(plan%sc)))
    plan%sc=0
    do p=1,size(slots)
       r=block_catalog(request(1,p))%owner+1
       pos=plan%sd(r)+plan%sc(r)+1
       send_key(1:4,pos)=request(:,p)
       send_key(5,pos)=grid(node_domain(slots(p)))%level%elts(node_id(slots(p))+1)
       send_key(6,pos)=node_flags(slots(p))
       plan%slot(pos)=slots(p)
       plan%sc(r)=plan%sc(r)+1
    end do
    call MPI_Alltoallv(send_key,6*plan%sc,6*plan%sd,MPI_INTEGER, &
         recv_key,6*plan%rc,6*plan%rd,MPI_INTEGER,comm,ierr)
    call check(ierr,'native dependency keys')
    plan%key=recv_key(1:4,1:sum(plan%rc))
    plan%source_info=send_key(5:6,1:sum(plan%sc))
    plan%target_info=recv_key(5:6,1:sum(plan%rc))
    call compile_local_block_inverse_routes(plan%key,plan%address)
  end subroutine

  subroutine select_transfer(plan,base,component,level,install,outer,modified_only)
    ! Select changed block levels plus every alias destination (all levels).
    ! A gather restores authoritative interiors overwritten by workspace
    ! aliases even when their block kernel did not run. Fixed-coarse resets
    ! likewise remain visible to boundary consumers on every phase.
    type(Transfer_Plan), intent(out) :: plan
    type(Transfer_Plan), intent(in) :: base
    integer, intent(in) :: component,level
    logical, intent(in) :: install,outer
    logical, optional, intent(in) :: modified_only
    integer :: r,p,s,t
    allocate(plan%sc(n_process),plan%sd(n_process),plan%rc(n_process),plan%rd(n_process))
    plan%sc=0
    plan%rc=0
    do r=1,n_process
       do p=base%sd(r)+1,base%sd(r)+base%sc(r)
          if (selected(base%source_info(:,p))) plan%sc(r)=plan%sc(r)+1
       end do
       do p=base%rd(r)+1,base%rd(r)+base%rc(r)
          if (selected(base%target_info(:,p))) plan%rc(r)=plan%rc(r)+1
       end do
    end do
    call displacements(plan%sc,plan%sd)
    call displacements(plan%rc,plan%rd)
    allocate(plan%slot(sum(plan%sc)),plan%key(4,sum(plan%rc)))
    if (present(modified_only)) then
       allocate(plan%source_info(2,sum(plan%sc)),plan%target_info(2,sum(plan%rc)))
    end if
    s=0
    do p=1,sum(base%sc)
       if (.not. selected(base%source_info(:,p))) cycle
       s=s+1
       plan%slot(s)=base%slot(p)
       if (present(modified_only)) plan%source_info(:,s)=base%source_info(:,p)
    end do
    t=0
    do p=1,sum(base%rc)
       if (.not. selected(base%target_info(:,p))) cycle
       t=t+1
       plan%key(:,t)=base%key(:,p)
       if (present(modified_only)) plan%target_info(:,t)=base%target_info(:,p)
    end do
    call compile_local_block_inverse_routes(plan%key,plan%address)
  contains
    logical function selected(info)
      integer, intent(in) :: info(2)
      if (install) then
         ! Remote aliases only write their requested level. Local aliases
         ! still write all levels; an outer transaction also writes coarse
         ! aliases before producing and synchronizing the fine level.
         selected=info(1)==level .or. btest(info(2),component+2) .or. btest(info(2),2)
         if (outer) selected=selected .or. info(1)==level-1
         if (present(modified_only)) then
            ! No block vector kernel intervenes since the preceding vector
            ! sync (or initial seeding). Unwritten interior/boundary records
            ! already contain the gathered value; retain aliases and resets.
            if (modified_only) selected=selected .and. &
                 (btest(info(2),1) .or. btest(info(2),2) .or. btest(info(2),5))
         end if
      else
         selected=info(1)==level .or. btest(info(2),component-1) .or. btest(info(2),2)
      end if
    end function
  end subroutine

  subroutine combine_transfers(plan,a,b)
    ! The two destinations are disjoint storage families (patch/boundary).
    ! Preserve each family's peer order, with one exchange and completion.
    ! Both ends derive this concatenation locally: no new setup collective.
    type(Transfer_Plan), intent(out) :: plan
    type(Transfer_Plan), intent(in) :: a,b
    integer :: r,s,t,na,nb
    plan%sc=a%sc+b%sc
    plan%rc=a%rc+b%rc
    allocate(plan%sd(n_process),plan%rd(n_process))
    call displacements(plan%sc,plan%sd)
    call displacements(plan%rc,plan%rd)
    allocate(plan%slot(sum(plan%sc)),plan%key(4,sum(plan%rc)))
    do r=1,n_process
       s=plan%sd(r)
       na=a%sc(r)
       nb=b%sc(r)
       plan%slot(s+1:s+na)=a%slot(a%sd(r)+1:a%sd(r)+na)
       plan%slot(s+na+1:s+na+nb)=b%slot(b%sd(r)+1:b%sd(r)+nb)
       t=plan%rd(r)
       na=a%rc(r)
       nb=b%rc(r)
       plan%key(:,t+1:t+na)=a%key(:,a%rd(r)+1:a%rd(r)+na)
       plan%key(:,t+na+1:t+na+nb)=b%key(:,b%rd(r)+1:b%rd(r)+nb)
    end do
    call compile_local_block_inverse_routes(plan%key,plan%address)
  end subroutine

  subroutine combine_outer_edges(plan,interior,boundary)
    ! A single peer exchange in edge units. Interior writes use the union of
    ! actual outer targets and signed alias destinations in both storage
    ! families. Scaffold resets retain all edges. Peer/family order is unchanged.
    type(Transfer_Plan), intent(out) :: plan
    type(Transfer_Plan), intent(in) :: interior,boundary
    integer :: r,p,e,s,t,n,pass
    integer, allocatable :: edges(:),address(:,:)
    plan%edgewise=.true.
    allocate(plan%sc(n_process),plan%sd(n_process),plan%rc(n_process),plan%rd(n_process))
    plan%sc=0
    plan%rc=0
    do pass=1,2
       s=0
       t=0
       do r=1,n_process
          do p=interior%sd(r)+1,interior%sd(r)+interior%sc(r)
             do e=1,EDGE
                if (.not.(btest(interior%source_info(2,p),5+e).or.btest(interior%source_info(2,p),2))) cycle
                s=s+1
                if (pass==1) plan%sc(r)=plan%sc(r)+1
                if (pass==2) plan%slot(s)=EDGE*(interior%slot(p)-1)+e
             end do
          end do
          do p=boundary%sd(r)+1,boundary%sd(r)+boundary%sc(r)
             do e=1,EDGE
                if (.not.(btest(boundary%source_info(2,p),5+e).or.btest(boundary%source_info(2,p),2))) cycle
                s=s+1
                if (pass==1) plan%sc(r)=plan%sc(r)+1
                if (pass==2) plan%slot(s)=EDGE*(boundary%slot(p)-1)+e
             end do
          end do
          do p=interior%rd(r)+1,interior%rd(r)+interior%rc(r)
             do e=1,EDGE
                if (.not.(btest(interior%target_info(2,p),5+e).or.btest(interior%target_info(2,p),2))) cycle
                t=t+1
                if (pass==1) plan%rc(r)=plan%rc(r)+1
                if (pass==2) then
                   plan%key(:,t)=interior%key(:,p)
                   edges(t)=e
                end if
             end do
          end do
          do p=boundary%rd(r)+1,boundary%rd(r)+boundary%rc(r)
             do e=1,EDGE
                if (.not.(btest(boundary%target_info(2,p),5+e).or.btest(boundary%target_info(2,p),2))) cycle
                t=t+1
                if (pass==1) plan%rc(r)=plan%rc(r)+1
                if (pass==2) then
                   plan%key(:,t)=boundary%key(:,p)
                   edges(t)=e
                end if
             end do
          end do
       end do
       if (pass==1) then
          call displacements(plan%sc,plan%sd)
          call displacements(plan%rc,plan%rd)
          allocate(plan%slot(s),plan%key(4,t),edges(t))
       end if
    end do
    call compile_local_block_inverse_routes(plan%key,address)
    n=size(address,2)
    allocate(plan%address(7,n))
    plan%address(1:6,:)=address
    plan%address(7,:)=edges
  end subroutine

  integer function alias_node(code,component) result(p)
    integer, intent(in) :: code,component
    p=abs(code)
    if (component==2) p=(p-1)/EDGE+1
  end function

  subroutine native_inverse_gather(component,family,reset_scaffold,changed_level)
    integer, intent(in) :: component,family
    logical, optional, intent(in) :: reset_scaffold
    integer, optional, intent(in) :: changed_level
    logical :: reset
    integer :: q,first,last
    reset=.true.
    if (present(reset_scaffold)) reset=reset_scaffold
    if (family==BLOCK_PAYLOAD_SOL .and. reset) then
       first=merge(1,nscalar+1,component==1)
       last=merge(nscalar,nscalar+EDGE,component==1)
       do q=1,size(scaffold_slot)
          value(first:last,:,scaffold_slot(q),family)=scaffold_value(first:last,:,q)
       end do
    end if
    if (resident203) then
       first=merge(1,nscalar+1,component==1);last=merge(nscalar,nscalar+EDGE,component==1)
       do q=1,size(restore203_slot)
          value(first:last,:,restore203_slot(q),family)=restore203(first:last,:,q)
       end do
       if (family/=BLOCK_PAYLOAD_SOL) call die('resident gather family')
       if (.not.present(changed_level)) call die('resident gather requires phase')
       if (oracle203) call check_gather203(level_interior(component,changed_level),component,family)
       if(ledger203) work203(5)=work203(5)+int(size(restore203_slot),int64)*(last-first+1)*inverse_levels
       if(ledger203) work203(1)=work203(1)+1_int64
       return
    end if
    if (present(changed_level)) then
       if (changed_level<level_start-2 .or. changed_level>level_end) call die('gather level range')
       call transfer_nodes(level_interior(component,changed_level),component,family,.false.)
    else
       call transfer_nodes(interior_plan,component,family,.false.)
    end if
  end subroutine

  subroutine native_inverse_scatter(component,family,interiors,changed_level)
    integer, intent(in) :: component,family
    logical, intent(in) :: interiors
    integer, optional, intent(in) :: changed_level
    if (resident203) then
       if (interiors) then
          if (.not.present(changed_level)) call die('resident outer install requires level')
          call remember_outer203(outer_install(changed_level))
       end if
       if(ledger203) work203(2)=work203(2)+1_int64
       if (.not.oracle203) return
    end if
    if (present(changed_level)) then
       if (changed_level<level_start-2 .or. changed_level>level_end) call die('scatter level range')
       if (interiors) then
          if (component/=2 .or. changed_level<level_start-1) call die('outer scatter phase')
          call transfer_nodes(outer_install(changed_level),component,family,.true.)
       else
          call transfer_nodes(level_boundary(component,changed_level),component,family,.true.)
       end if
    else
       if (interiors) call transfer_nodes(interior_plan,component,family,.true.)
       call transfer_nodes(boundary_plan,component,family,.true.)
    end if
  end subroutine

#ifdef STAGE202_INVERSE
  subroutine prepare_self_routes(plan)
    type(Transfer_Plan), intent(inout) :: plan
    integer :: first,last
    first=plan%sd(rank+1)+1
    last=plan%sd(rank+1)+plan%sc(rank+1)
    if (plan%sc(rank+1)/=plan%rc(rank+1)) call die('native self route extent')
    plan%self_slot=plan%slot(first:last)
    if (plan%edgewise) plan%self_slot=(plan%self_slot-1)/EDGE+1
  end subroutine prepare_self_routes
#endif

  subroutine transfer_nodes(plan,component,family,install)
    type(Transfer_Plan), intent(in) :: plan
    integer, intent(in) :: component,family
    logical, intent(in) :: install
    integer :: p,q,k,v,nv,first,n,phase
#ifdef STAGE202_INVERSE
    integer :: node,source_first,self_first,self_last,local_first
#endif
    real(dp) :: started
    first=merge(1,nscalar+1,component==1)
    nv=merge(nscalar,EDGE,component==1)
    if (plan%edgewise) then
       if (.not.install .or. component/=2) call die('edgewise transfer outside vector install')
       nv=1
    end if
    n=nv*inverse_levels
#ifdef STAGE202_INVERSE
    self_first=plan%sd(rank+1)+1
    self_last=plan%sd(rank+1)+plan%sc(rank+1)
    local_first=plan%rd(rank+1)+1
    ! Request order on the self peer pairs these workspace slots with the
    ! corresponding block routes. Edge installs encode component in slot.
#endif
    started=0.0_dp
    if (timing) started=MPI_Wtime()
    if (install) then
       call detail_enter(DP_INVERSE_PACK)
       do p=1,sum(plan%sc)
#ifdef STAGE202_INVERSE
          if (p>=self_first.and.p<=self_last) cycle
#endif
          q=(p-1)*n
#ifdef STAGE202_INVERSE
          node=plan%slot(p)
          source_first=first
          if (plan%edgewise) then
             source_first=nscalar+mod(node-1,EDGE)+1
             node=(node-1)/EDGE+1
          end if
          do k=1,inverse_levels
             do v=1,nv
                q=q+1
                transfer_send(q)=value(source_first+v-1,k,node,family)
             end do
          end do
#else
          do k=1,inverse_levels
             do v=1,nv
                q=q+1
                if (plan%edgewise) then
                   transfer_send(q)=value(nscalar+mod(plan%slot(p)-1,EDGE)+1,k,(plan%slot(p)-1)/EDGE+1,family)
                else
                   transfer_send(q)=value(first+v-1,k,plan%slot(p),family)
                end if
             end do
          end do
#endif
       end do
       call detail_leave(DP_INVERSE_PACK)

#ifdef STAGE202_INVERSE
       call exchange_values(transfer_send,plan%sc,plan%sd,transfer_recv,plan%rc,plan%rd,n,19071,.false.)
#else
       call exchange_values(transfer_send,plan%sc,plan%sd,transfer_recv,plan%rc,plan%rd,n,19071)
#endif
       call detail_enter(DP_INVERSE_INSTALL)
       if (plan%edgewise) then
#ifdef STAGE202_INVERSE
          call transfer_local_block_inverse_routes(plan%address(1:6,:),family,component,.true.,transfer_recv, &
               plan%address(7,:),value(:,:,:,family),local_first,plan%self_slot, &
            field_first=inverse_field_first,field_count=inverse_levels)
#else
       call transfer_local_block_inverse_routes(plan%address(1:6,:),family,component,.true.,transfer_recv,plan%address(7,:), &
            field_first=inverse_field_first,field_count=inverse_levels)
#endif
       else
#ifdef STAGE202_INVERSE
          call transfer_local_block_inverse_routes(plan%address,family,component,.true.,transfer_recv, &
               workspace=value(:,:,:,family),local_first=local_first,local_slots=plan%self_slot, &
            field_first=inverse_field_first,field_count=inverse_levels)
#else
       call transfer_local_block_inverse_routes(plan%address,family,component,.true.,transfer_recv, &
            field_first=inverse_field_first,field_count=inverse_levels)
#endif
       end if
       call detail_leave(DP_INVERSE_INSTALL)
    else
       call detail_enter(DP_INVERSE_PACK)

#ifdef STAGE202_INVERSE
       call transfer_local_block_inverse_routes(plan%address,family,component,.false.,transfer_recv, &
            workspace=value(:,:,:,family),local_first=local_first,local_slots=plan%self_slot, &
            field_first=inverse_field_first,field_count=inverse_levels)
#else
       call transfer_local_block_inverse_routes(plan%address,family,component,.false.,transfer_recv, &
            field_first=inverse_field_first,field_count=inverse_levels)
#endif
       call detail_leave(DP_INVERSE_PACK)

#ifdef STAGE202_INVERSE
       call exchange_values(transfer_recv,plan%rc,plan%rd,transfer_send,plan%sc,plan%sd,n,19072,.false.)
#else
       call exchange_values(transfer_recv,plan%rc,plan%rd,transfer_send,plan%sc,plan%sd,n,19072)
#endif
       call detail_enter(DP_INVERSE_INSTALL)
       do p=1,sum(plan%sc)
#ifdef STAGE202_INVERSE
          if (p>=self_first.and.p<=self_last) cycle
#endif
          q=(p-1)*n
          do k=1,inverse_levels
             do v=1,nv
                q=q+1
                value(first+v-1,k,plan%slot(p),family)=transfer_send(q)
             end do
          end do
       end do
       call detail_leave(DP_INVERSE_INSTALL)
    end if
    if (timing) then
       phase=merge(3,1,install)
       native_inverse_seconds(phase)=native_inverse_seconds(phase)+MPI_Wtime()-started
       native_inverse_calls(phase)=native_inverse_calls(phase)+1_int64
       if (install) then
          call record_traffic(phase,plan%sc,n)
       else
          call record_traffic(phase,plan%rc,n)
       end if
    end if
  end subroutine

  subroutine exchange_values(send,sc,sd,recv,rc,rd,n,tag,copy_self)
    real(dp), intent(in), asynchronous :: send(:)
    real(dp), intent(out), asynchronous :: recv(:)
    integer, intent(in) :: sc(:),sd(:),rc(:),rd(:),n,tag
    logical, optional, intent(in) :: copy_self
    integer :: nrequest
    call begin_exchange_values(send,sc,sd,recv,rc,rd,n,tag,nrequest,copy_self)
    call finish_exchange_values(nrequest)
  end subroutine

  subroutine begin_exchange_values(send,sc,sd,recv,rc,rd,n,tag,nrequest,copy_self)
    ! Persistent sparse peer schedules, with no per-phase global collective.
    real(dp), intent(in), asynchronous :: send(:)
    real(dp), intent(out), asynchronous :: recv(:)
    integer, intent(in) :: sc(:),sd(:),rc(:),rd(:),n,tag
    integer, intent(out) :: nrequest
    logical, optional, intent(in) :: copy_self
    logical :: self_copy
    integer :: r,ierr
    call detail_enter(DP_INVERSE_POST)
    nrequest=0
    do r=1,n_process
       if (r==rank+1 .or. rc(r)==0) cycle
       nrequest=nrequest+1
       call MPI_Irecv(recv(n*rd(r)+1:n*(rd(r)+rc(r))),n*rc(r),MPI_DOUBLE_PRECISION, &
            r-1,tag,comm,requests(nrequest),ierr)
       call check(ierr,'native dependency receive')
    end do
    do r=1,n_process
       if (r==rank+1 .or. sc(r)==0) cycle
       nrequest=nrequest+1
       call MPI_Isend(send(n*sd(r)+1:n*(sd(r)+sc(r))),n*sc(r),MPI_DOUBLE_PRECISION, &
            r-1,tag,comm,requests(nrequest),ierr)
       call check(ierr,'native dependency send')
    end do
    r=rank+1
    if (sc(r)/=rc(r)) call die('native self route extent')
    self_copy=.true.
    if (present(copy_self)) self_copy=copy_self
    if (self_copy) then
    recv(n*rd(r)+1:n*(rd(r)+rc(r)))=send(n*sd(r)+1:n*(sd(r)+sc(r)))
    end if
    call detail_leave(DP_INVERSE_POST)

  end subroutine

  subroutine finish_exchange_values(nrequest)
    integer, intent(in) :: nrequest
    integer :: ierr
    call detail_enter(DP_INVERSE_WAIT)
    if (nrequest>0) then
       call MPI_Waitall(nrequest,requests(1:nrequest),MPI_STATUSES_IGNORE,ierr)
       call check(ierr,'native dependency completion')
    end if
    call detail_leave(DP_INVERSE_WAIT)

  end subroutine

  subroutine record_traffic(phase,counts,nvalue)
    integer, intent(in) :: phase,counts(:),nvalue
    integer :: r
    do r=1,n_process
       if (r==rank+1 .or. counts(r)==0) cycle
       native_inverse_messages(phase)=native_inverse_messages(phase)+1_int64
       native_inverse_bytes(phase)=native_inverse_bytes(phase)+8_int64*nvalue*counts(r)
    end do
  end subroutine

  subroutine compile_aliases
    integer :: kind,pos,l,r,ds,dd,g,q,id,dest,n,s,t
    do kind=1,2
       pos=merge(AT_NODE,AT_EDGE,kind==1)
       do l=level_start-1,level_end
          associate(a=>aliases(kind,l))
            allocate(a%sc(n_process),a%sd(n_process),a%rc(n_process),a%rd(n_process))
            a%sc=0
            a%rc=0
            do r=1,n_process
               if (r==rank+1) cycle
               do ds=1,size(grid)
                  do dd=1,n_domain(r)
                     g=glo_id(r,dd)+1
                     do q=1,grid(ds)%pack(pos,g)%length
                        id=grid(ds)%pack(pos,g)%elts(q)
                        if (grid(ds)%level%elts(id/merge(1,EDGE,kind==1)+1)==l) a%sc(r)=a%sc(r)+1
                     end do
                  end do
               end do
               do ds=1,n_domain(r)
                  g=glo_id(r,ds)+1
                  do dd=1,size(grid)
                     do q=1,grid(dd)%unpk(pos,g)%length
                        id=abs(grid(dd)%unpk(pos,g)%elts(q))
                        if (grid(dd)%level%elts(id/merge(1,EDGE,kind==1)+1)==l) a%rc(r)=a%rc(r)+1
                     end do
                  end do
               end do
            end do
            call displacements(a%sc,a%sd)
            call displacements(a%rc,a%rd)
            allocate(a%source(sum(a%sc)),a%dest(sum(a%rc)))
            s=0
            t=0
            do r=1,n_process
               if (r==rank+1) cycle
               do ds=1,size(grid)
                  do dd=1,n_domain(r)
                     g=glo_id(r,dd)+1
                     do q=1,grid(ds)%pack(pos,g)%length
                        id=grid(ds)%pack(pos,g)%elts(q)
                        if (grid(ds)%level%elts(id/merge(1,EDGE,kind==1)+1)/=l) cycle
                        s=s+1
                        a%source(s)=route_slot(ds,id,kind)
                     end do
                  end do
               end do
               do ds=1,n_domain(r)
                  g=glo_id(r,ds)+1
                  do dd=1,size(grid)
                     do q=1,grid(dd)%unpk(pos,g)%length
                        id=grid(dd)%unpk(pos,g)%elts(q)
                        if (grid(dd)%level%elts(abs(id)/merge(1,EDGE,kind==1)+1)/=l) cycle
                        t=t+1
                        a%dest(t)=route_slot(dd,id,kind)
                     end do
                  end do
               end do
            end do
          end associate
       end do
       n=0
       do ds=1,size(grid)
          do dd=1,size(grid)
             n=n+grid(ds)%pack(pos,glo_id(rank+1,dd)+1)%length
          end do
       end do
       allocate(local_alias(kind)%source(n),local_alias(kind)%dest(n))
       n=0
       do ds=1,size(grid)
          do dd=1,size(grid)
             g=glo_id(rank+1,dd)+1
             dest=glo_id(rank+1,ds)+1
             if (grid(ds)%pack(pos,g)%length /= grid(dd)%unpk(pos,dest)%length) call die('local alias extent')
             do q=1,grid(ds)%pack(pos,g)%length
                n=n+1
                local_alias(kind)%source(n)=route_slot(ds,grid(ds)%pack(pos,g)%elts(q),kind)
                local_alias(kind)%dest(n)=route_slot(dd,grid(dd)%unpk(pos,dest)%elts(q),kind)
             end do
          end do
       end do
    end do
  end subroutine

  subroutine combine_aliases(a,component,first)
    ! Pack all requested levels into one peer message, but retain the old
    ! level-major receive installation order, including repeated targets.
    type(Alias_Plan), intent(out) :: a
    integer, intent(in) :: component,first
    integer :: l,r,p,s,t,q
    integer :: cursor(n_process)
    allocate(a%sc(n_process),a%sd(n_process),a%rc(n_process),a%rd(n_process))
    a%sc=0
    a%rc=0
    do l=first,level_end
       a%sc=a%sc+aliases(component,l)%sc
       a%rc=a%rc+aliases(component,l)%rc
    end do
    call displacements(a%sc,a%sd)
    call displacements(a%rc,a%rd)
    allocate(a%source(sum(a%sc)),a%dest(sum(a%rc)),a%install_order(sum(a%rc)))
    s=0
    t=0
    do r=1,n_process
       do l=first,level_end
          associate(b=>aliases(component,l))
            do p=b%sd(r)+1,b%sd(r)+b%sc(r)
               s=s+1
               a%source(s)=b%source(p)
            end do
            do p=b%rd(r)+1,b%rd(r)+b%rc(r)
               t=t+1
               a%dest(t)=b%dest(p)
            end do
          end associate
       end do
    end do
    cursor=a%rd
    q=0
    do l=first,level_end
       do r=1,n_process
          do p=1,aliases(component,l)%rc(r)
             q=q+1
             cursor(r)=cursor(r)+1
             a%install_order(q)=cursor(r)
          end do
       end do
    end do
  end subroutine

  subroutine native_inverse_boundary(component,family,first,last)
    integer, intent(in) :: component,family,first,last
    real(dp) :: started
    started=0.0_dp
    if (timing) started=MPI_Wtime()
    if (first==last) then
       call apply_aliases(aliases(component,first),component,family)
    else if (last==level_end .and. first==level_start-1) then
       call apply_aliases(range_alias(component,1),component,family)
    else if (last==level_end .and. first==level_start) then
       call apply_aliases(range_alias(component,2),component,family)
    else
       call die('uncompiled alias level range')
    end if
    if (timing) then
       native_inverse_seconds(2)=native_inverse_seconds(2)+MPI_Wtime()-started
       native_inverse_calls(2)=native_inverse_calls(2)+1_int64
    end if
  end subroutine

  subroutine apply_aliases(a,component,family)
    type(Alias_Plan), intent(inout) :: a
    integer, intent(in) :: component,family
    integer :: p,k,v,nv,n,q,t,nrequest
#ifdef STAGE202_INVERSE
    integer :: sn,dn,sf,df
    logical :: negative
#endif
    nv=merge(nscalar,1,component==1)
    n=nv*inverse_levels
    ! Remote sends are snapshots BEFORE the ordered, all-level local copies.
    call detail_enter(DP_INVERSE_PACK)
    q=0
    do p=1,size(a%source)
#ifdef STAGE202_INVERSE
       call alias_address(a%source(p),component,sn,sf)
       do k=1,inverse_levels
          do v=1,nv
             q=q+1
             alias_send(q)=value(sf+v-1,k,sn,family)
          end do
       end do
#else
       do k=1,inverse_levels
          do v=1,nv
             q=q+1
             alias_send(q)=alias_value(a%source(p),component,v,k,family)
          end do
       end do
#endif
    end do
    ! Sends have already snapshotted every remote source. Local copies use
    ! only value(), so communication can progress without changing ordering.
    call detail_leave(DP_INVERSE_PACK)
    call begin_exchange_values(alias_send,a%sc,a%sd,alias_recv,a%rc,a%rd,n,19073,nrequest)
    call detail_enter(DP_INVERSE_LOCAL)
    do p=1,size(local_alias(component)%source)
#ifdef STAGE202_INVERSE
       call alias_address(local_alias(component)%source(p),component,sn,sf)
       call alias_address(local_alias(component)%dest(p),component,dn,df)
       negative=component==2.and.local_alias(component)%dest(p)<0
       ! Keep p,k,v ordering: local alias chains can consume earlier writes.
       if (negative) then
          do k=1,inverse_levels
             do v=1,nv
                value(df+v-1,k,dn,family)=-value(sf+v-1,k,sn,family)
             end do
          end do
       else
          do k=1,inverse_levels
             do v=1,nv
                value(df+v-1,k,dn,family)=value(sf+v-1,k,sn,family)
             end do
          end do
       end if
#else
       do k=1,inverse_levels
          do v=1,nv
             call set_alias(local_alias(component)%dest(p),component,v,k,family, &
                  alias_value(local_alias(component)%source(p),component,v,k,family))
          end do
       end do
#endif
    end do
    call detail_leave(DP_INVERSE_LOCAL)
    call finish_exchange_values(nrequest)
    if (timing) call record_traffic(2,a%sc,n)
    call detail_enter(DP_INVERSE_INSTALL)
    do p=1,size(a%dest)
       t=p
       if (allocated(a%install_order)) t=a%install_order(p)
       q=(t-1)*n
#ifdef STAGE202_INVERSE
       call alias_address(a%dest(t),component,dn,df)
       if (component==2.and.a%dest(t)<0) then
          do k=1,inverse_levels
             do v=1,nv
                q=q+1
                value(df+v-1,k,dn,family)=-alias_recv(q)
             end do
          end do
       else
          do k=1,inverse_levels
             do v=1,nv
                q=q+1
                value(df+v-1,k,dn,family)=alias_recv(q)
             end do
          end do
       end if
#else
       do k=1,inverse_levels
          do v=1,nv
             q=q+1
             call set_alias(a%dest(t),component,v,k,family,alias_recv(q))
          end do
       end do
#endif
    end do
    call detail_leave(DP_INVERSE_INSTALL)
  end subroutine

#ifdef STAGE202_INVERSE
  subroutine alias_address(code,component,node,first)
    integer, intent(in) :: code,component
    integer, intent(out) :: node,first
    node=code
    first=1
    if (component==2) then
       node=(abs(code)-1)/EDGE+1
       first=nscalar+mod(abs(code)-1,EDGE)+1
    end if
  end subroutine alias_address
#endif

  real(dp) function alias_value(code,component,v,k,family) result(x)
    integer, intent(in) :: code,component,v,k,family
    if (component==1) then
       x=value(v,k,code,family)
    else
       x=value(nscalar+mod(code-1,EDGE)+1,k,(code-1)/EDGE+1,family)
    end if
  end function

  subroutine set_alias(code,component,v,k,family,x)
    integer, intent(in) :: code,component,v,k,family
    real(dp), intent(in) :: x
    integer :: q
    if (component==1) then
       value(v,k,code,family)=x
    else
       q=abs(code)
       value(nscalar+mod(q-1,EDGE)+1,k,(q-1)/EDGE+1,family)=x
       if (code < 0) value(nscalar+mod(q-1,EDGE)+1,k,(q-1)/EDGE+1,family)=-x
    end if
  end subroutine

  subroutine compare_native_inverse(scaling,component,payload_family,consumer_aliases,exact)
    type(Float_Field), intent(in) :: scaling(:,:)
    logical,intent(in),optional::consumer_aliases,exact
    logical::consumer,strict
    integer, intent(in) :: component
    integer, optional, intent(in) :: payload_family
    integer :: p,d,id,k,v,first,last,var,address,family
    real(dp) :: reference,native,allowed
    consumer=.false.;strict=.false.
    if(present(consumer_aliases)) consumer=consumer_aliases
    if(present(exact)) strict=exact
    if(size(scaling,1)/=N_VARIABLE.or.size(scaling,2)/=inverse_levels) call die('inverse comparison shape')
    first=merge(1,nscalar+1,component==1)
    last=merge(nscalar,nscalar+EDGE,component==1)
    family=BLOCK_PAYLOAD_SOL
    if (present(payload_family)) family=payload_family
    do p=1,nnode
       ! Additional registered nodes only supply fixed scaffold staging; they
       ! were not published to the independent provisional reference field.
       if(block_rk_enabled.and.inverse_field_first==1.and.inverse_levels==zlevels.and.p>original_inverse_nodes) cycle
       d=node_domain(p)
       id=node_id(p)
       do k=1,inverse_levels
          do v=first,last
             if (component==1) then
                var=scalars(1)+v-1
                address=id+1
             else
                var=S_VELO
                address=EDGE*id+v-nscalar
             end if
             native=value(v,k,p,family)
             if(consumer.and.nodes(d)%key(1,id+1)==0) &
                  native=consumer_alias_value(v,k,consumer_alias_index(p))
             reference=scaling(var,k)%data(d)%elts(address)
             if (transfer_bits(native)==transfer_bits(reference)) cycle
             if (.not. ieee_is_finite(native) .or. .not. ieee_is_finite(reference)) call die('nonfinite oracle difference')
             allowed=128.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(native),abs(reference))
             if (.not.strict.and.abs(native-reference) <= allowed) cycle
             write(error_unit,'(a,5i8,3es25.16)') 'Inverse phase mismatch d,node,k,var,component: ', &
                  d,id,k,var,v,native,reference,native-reference
             write(error_unit,'(a,5i8)') 'level and native key: ',grid(d)%level%elts(id+1),nodes(d)%key(:,id+1)
             call die('phase oracle mismatch')
          end do
       end do
    end do
  end subroutine

  integer(int64) function transfer_bits(x) result(bits)
    real(dp), intent(in) :: x
    bits=transfer(x,bits)
  end function

  subroutine physics_consumer_alias(d,id,operand,mean,soil)
    ! The saved final-inverse consumer version is distinct from the current
    ! initial/provisional tendency epoch and from canonical block interiors.
    integer,intent(in)::d,id
    real(dp),intent(out)::operand(nscalar+EDGE,zlevels)
    real(dp),intent(out),optional::mean(nscalar,zlevels),soil(1-zmin)
    integer::key(4),p,q
    if(inverse_field_first/=1.or.inverse_levels/=zlevels) call die('physics accessor received non-atmospheric epoch')
    if(.not.physics_consumer_ready) call die('physics consumer alias epoch is not ready')
    call inverse_operand_key(d,id,key)
    if(key(1)/=0) call die('physics alias accessor attempted canonical interior')
    p=nodes(d)%slot(id+1)
    if(p<1.or.p>original_inverse_nodes) call die('physics alias is outside the consumer closure')
    q=consumer_alias_index(p)
    if(q<1) call die('physics alias index absent')
    if((present(mean).or.present(soil)).and..not.physics_reference_ready) &
         call die('physics reference alias epoch is not ready')
    operand=consumer_alias_value(:,:,q)
    if(present(mean)) mean=physics_alias_mean(:,:,q)
    if(present(soil)) soil=physics_alias_soil(S_TEMP-scalars(1)+1,:,q)
    if(.not.all(ieee_is_finite(operand))) call die('nonfinite saved physics alias')
  end subroutine

  subroutine install_physics_consumer_alias(d,id,result,write_velocity)
    integer,intent(in)::d,id
    real(dp),intent(in)::result(4*zlevels+1-zmin)
    logical,intent(in)::write_velocity
    real(dp)::entry(nscalar+EDGE,zlevels)
    integer::q,k
    call physics_consumer_alias(d,id,entry)
    q=consumer_alias_index(nodes(d)%slot(id+1))
    consumer_alias_value(S_TEMP-scalars(1)+1,:,q)=result(3*zlevels+1:4*zlevels)
    if(write_velocity) then
       do k=1,zlevels
          consumer_alias_value(nscalar+1:nscalar+EDGE,k,q)=result(3*(k-1)+1:3*k)
       end do
    end if
    physics_alias_soil(S_TEMP-scalars(1)+1,:,q)=result(4*zlevels+1:)
  end subroutine install_physics_consumer_alias

  subroutine export_native_soil_aliases(address,operand,generation)
    integer,allocatable,intent(out)::address(:,:)
    real(dp),allocatable,intent(out)::operand(:,:,:)
    integer(int64),intent(out)::generation
    integer::q,p
    if(.not.physics_results_ready.or..not.physics_reference_ready) call die('post-physics soil epoch absent')
    if(inverse_field_first/=1.or.inverse_levels/=zlevels) call die('soil export requires completed atmosphere')
    allocate(address(2,size(consumer_alias_slot)),operand(nscalar+EDGE,1-zmin,size(consumer_alias_slot)))
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q);address(:,q)=[node_domain(p),node_id(p)]
       operand(:,:,q)=physics_alias_soil(:,1-zmin:1:-1,q)
    end do
    generation=plan_generation
  end subroutine export_native_soil_aliases

  subroutine finish_physics_consumer_epoch
    if(.not.physics_consumer_ready) call die('physics consumer epoch already consumed')
    physics_consumer_ready=.false.
#ifdef PHYSICS_NATIVE_STATE
    physics_results_ready=.true.
#endif
  end subroutine

  subroutine capture_native_inverse_consumer_aliases(component,level)
    integer, intent(in) :: component,level
    integer :: q,p,d,id,first,last
    first=merge(1,nscalar+1,component==1)
    last=merge(nscalar,nscalar+EDGE,component==1)
    do q=1,size(consumer_alias_slot)
       p=consumer_alias_slot(q)
       d=node_domain(p)
       id=node_id(p)
       if (grid(d)%level%elts(id+1)/=level) cycle
       consumer_alias_value(first:last,:,q)=value(first:last,:,p,BLOCK_PAYLOAD_SOL)
    end do
  end subroutine

  subroutine publish_native_inverse(scaling,payload_family,consumer_aliases,all_nodes)
    ! One compatibility publication after the complete native transaction.
    ! Interior installation remains the standard final-owner publication.
    type(Float_Field), intent(inout) :: scaling(:,:)
    integer, optional, intent(in) :: payload_family
    logical, optional, intent(in) :: consumer_aliases,all_nodes
    integer :: p,d,id,k,v,family,q
    logical :: consumer,publish_all
    publish_all=.false.
    if(present(all_nodes)) publish_all=all_nodes
    if(size(scaling,1)/=N_VARIABLE.or.size(scaling,2)/=inverse_levels) call die('inverse publication shape')
    family=BLOCK_PAYLOAD_SOL
    if (present(payload_family)) family=payload_family
    consumer=.false.
    if (present(consumer_aliases)) consumer=consumer_aliases
    if (consumer.and.family/=BLOCK_PAYLOAD_SOL) call die('consumer aliases require solution')
    if (consumer) then
       physics_consumer_ready=inverse_field_first==1.and.inverse_levels==zlevels
       do q=1,size(consumer_alias_slot)
          p=consumer_alias_slot(q)
          if(block_rk_enabled.and.inverse_field_first==1.and.inverse_levels==zlevels.and.p>original_inverse_nodes) cycle
          d=node_domain(p)
          id=node_id(p)
          do k=1,inverse_levels
             do v=1,nscalar
                scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)=consumer_alias_value(v,k,q)
             end do
             do v=1,EDGE
                scaling(S_VELO,k)%data(d)%elts(EDGE*id+v)=consumer_alias_value(nscalar+v,k,q)
             end do
          end do
       end do
       return
    end if
    do p=1,nnode
       if(block_rk_enabled.and.inverse_field_first==1.and.inverse_levels==zlevels.and.p>original_inverse_nodes) cycle
       d=node_domain(p)
       id=node_id(p)
       if (nodes(d)%key(1,id+1)>0.and..not.publish_all) cycle
       do k=1,inverse_levels
          do v=1,nscalar
             scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)=value(v,k,p,family)
          end do
          do v=1,EDGE
             scaling(S_VELO,k)%data(d)%elts(EDGE*id+v)=value(nscalar+v,k,p,family)
          end do
       end do
    end do
  end subroutine

  subroutine compile_outer
    integer :: l,d,p,pi,c,child,i,j,ic,jc,ip,jp,e,corner,q,a,b,t,mode
    integer :: op(N_BDRY+1),oc(N_BDRY+1),dp0(2,N_BDRY+1),dc(2,N_BDRY+1)
    type(Outer_Operation) :: item
    do l=level_start-1,level_end-1
       operation_first(l)=noperation+1
       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             if (.not. any(grid(d)%patch%elts(p+1)%children>0)) cycle
             coverage(1,l)=coverage(1,l)+int(inverse_levels,int64)
             call get_offs_Domain(grid(d),p,op,dp0)
             do c=1,N_CHDRN
                child=grid(d)%patch%elts(p+1)%children(c)
                if (child==0) cycle
                call get_offs_Domain(grid(d),child,oc,dc)
                do j=1,PATCH_SIZE/2+1
                   jc=2*(j-1)
                   jp=j-1+chd_offs(2,c)
                   do i=1,PATCH_SIZE/2+1
                      ic=2*(i-1)
                      ip=i-1+chd_offs(1,c)
                      coverage(2,l)=coverage(2,l)+int(inverse_levels,int64)
                      do e=RT,UP
                         item=Outer_Operation(level=l)
                         a=idx(ic+end_pt(1,1,e+1),jc+end_pt(2,1,e+1),oc,dc)
                         b=idx(ic+end_pt(1,2,e+1),jc+end_pt(2,2,e+1),oc,dc)
                         item%target=[edge_slot(d,EDGE*idx(ip,jp,op,dp0)+e), &
                              edge_slot(d,EDGE*a+e),edge_slot(d,EDGE*b+e)]
                         item%weight=base_weight+[(grid(d)%I_u_wgt%elts(b+1)%enc(q),q=1,9)]
                         item%source(1)=item%target(1)
                         item%source(2)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+3)
                         item%source(3)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+4)
                         item%source(4)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+6)
                         item%source(5)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+1)
                         item%source(6)=side(ip+opp_no(1,1,e+1),jp+opp_no(2,1,e+1),hex_s_offs(e+1)+2)
                         item%source(7)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+3)
                         item%source(8)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+4)
                         item%source(9)=side(ip+opp_no(1,1,e+1),jp+opp_no(2,1,e+1),hex_s_offs(e+1)+5)
                         item%source(10)=side(ip+opp_no(1,2,e+1),jp+opp_no(2,2,e+1),hex_s_offs(e+1)+5)
                         item%source(11)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+6)
                         item%source(12)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+1)
                         item%source(13)=side(ip+opp_no(1,2,e+1),jp+opp_no(2,2,e+1),hex_s_offs(e+1)+2)
                         call append_operation
                         coverage(3,l)=coverage(3,l)+2_int64*inverse_levels
                      end do
                   end do
                end do
             end do
          end do
          ! Preserve the original per-Domain correction order, after all
          ! regular plus-side writes on that Domain and before the next one.
          do corner=NORTHEAST,NORTHWEST
             if (.not. grid(d)%penta(corner)) cycle
             p=1
             do while (p>0)
                child=grid(d)%patch%elts(p+1)%children(corner-4)
                if (grid(d)%patch%elts(p+1)%level<l) then
                   p=child
                   cycle
                end if
                if (grid(d)%patch%elts(p+1)%level>l) exit
                if (child>0) then
                   call get_offs_Domain(grid(d),p,op,dp0)
                   call get_offs_Domain(grid(d),child,oc,dc)
                   if (corner==IMINUSJPLUS) then
                      item=Outer_Operation(level=l,kind=1)
                      a=idx(0,LAST-1,oc,dc)
                      b=idx(0,LAST,oc,dc)
                      item%target(2:3)=[edge_slot(d,EDGE*a+UP),edge_slot(d,EDGE*b+UP)]
                      item%source(1:2)=[edge_slot(d,EDGE*idx(0,PATCH_SIZE,op,dp0)+UP), &
                           edge_slot(d,EDGE*idx(-1,PATCH_SIZE,op,dp0)+RT)]
                      item%weight(1)=base_weight(8)+grid(d)%I_u_wgt%elts(b+1)%enc(8)
                      call append_operation
                   else if (corner==IPLUSJMINUS) then
                      item=Outer_Operation(level=l,kind=1)
                      a=idx(LAST-1,0,oc,dc)
                      b=idx(LAST,0,oc,dc)
                      item%target(2:3)=[edge_slot(d,EDGE*a+RT),edge_slot(d,EDGE*b+RT)]
                      item%source(1:2)=[edge_slot(d,EDGE*idx(PATCH_SIZE,0,op,dp0)+RT), &
                           edge_slot(d,EDGE*idx(PATCH_SIZE,-1,op,dp0)+UP)]
                      item%weight(1)=-(base_weight(7)+grid(d)%I_u_wgt%elts(b+1)%enc(7))
                      call append_operation
                   else if (corner==IJMINUS) then
                      do mode=2,3
                         item=Outer_Operation(level=l,kind=mode)
                         e=merge(UP,RT,mode==2)
                         a=idx(0,0,oc,dc)
                         b=idx(merge(0,1,mode==2),merge(1,0,mode==2),oc,dc)
                         item%target(2:3)=[edge_slot(d,EDGE*a+e),edge_slot(d,EDGE*b+e)]
                         t=idx(end_pt(1,2,e+1),end_pt(2,2,e+1),oc,dc)
                         q=merge(9,6,mode==2)
                         item%weight(1)=base_weight(q)+grid(d)%I_u_wgt%elts(t+1)%enc(q)
                         if (mode==2) then
                            item%source(1:2)=[edge_slot(d,EDGE*idx(0,-1,op,dp0)+UP), &
                                 edge_slot(d,EDGE*idx(-1,-1,op,dp0)+RT)]
                            item%source(3)=side(end_pt(1,1,UP+1),end_pt(2,1,UP+1),hex_s_offs(UP+1)+1)
                            item%source(4)=side(opp_no(1,2,UP+1),opp_no(2,2,UP+1),hex_s_offs(UP+1)+2)
                         else
                            item%source(1:2)=[edge_slot(d,EDGE*idx(-1,-1,op,dp0)+RT), &
                                 edge_slot(d,EDGE*idx(-1,0,op,dp0)+RT)]
                            item%source(3)=side(opp_no(1,1,RT+1),opp_no(2,1,RT+1),hex_s_offs(RT+1)+2)
                            item%source(4)=side(end_pt(1,1,RT+1),end_pt(2,1,RT+1),hex_s_offs(RT+1)+3)
                         end if
                         call append_operation
                      end do
                   end if
                end if
                p=child
             end do
          end do
       end do
       operation_last(l)=noperation
    end do
  contains
    integer function side(i,j,s) result(address)
      integer, intent(in) :: i,j,s
      address=edge_slot(d,ed_idx(i,j,hex_sides(:,s),op,dp0))
    end function
    subroutine append_operation
      noperation=noperation+1
      if (noperation>size(operations)) call die('outer tape capacity')
      operations(noperation)=item
      if (item%kind/=0) coverage(4,l)=coverage(4,l)+2_int64*inverse_levels
    end subroutine
  end subroutine

  subroutine native_inverse_outer(level,count)
    integer, intent(in) :: level
    integer(int64), intent(out) :: count(4)
    integer :: p,k,q
    real(dp) :: s(13),x(9),correction,second,first
    call detail_enter(DP_OUTER_KERNEL)
    count=coverage(:,level)
    first=0.0_dp
    second=0.0_dp
    do k=1,inverse_levels
       do p=operation_first(level),operation_last(level)
          associate(op=>operations(p))
            do q=1,13
               if (op%source(q)==0) exit
               s(q)=alias_value(op%source(q),2,1,k,1)
            end do
            select case(op%kind)
            case(0)
               x(1:5)=s(1:5)
               x(6:9)=[s(6)-s(7),s(8)-s(9),s(10)-s(11),s(12)-s(13)]
               second=sum(op%weight*x)+alias_value(op%target(3),2,1,k,2)
               call set_alias(op%target(3),2,1,k,1,second)
               ! Read the parent after the second write, matching the exact
               ! in-place dependency order even for aliased geometry slots.
               first=2.0_dp*alias_value(op%target(1),2,1,k,1)-second
            case(1:3)
               select case(op%kind)
               case(1)
                  correction=op%weight(1)*(s(1)+s(2))
               case(2)
                  correction=op%weight(1)*((-s(1)-(-s(2)))-(s(3)-s(4)))
               case(3)
                  correction=op%weight(1)*(s(1)+s(2)-(s(3)-s(4)))
               end select
               first=alias_value(op%target(2),2,1,k,1)-correction
               call set_alias(op%target(2),2,1,k,1,first)
               second=alias_value(op%target(3),2,1,k,1)+correction
               call set_alias(op%target(3),2,1,k,1,second)
            case default
               call die('invalid outer operation')
            end select
            if (.not. ieee_is_finite(first) .or. .not. ieee_is_finite(second)) call die('nonfinite outer edge')
            call set_alias(op%target(2),2,1,k,1,first)
          end associate
       end do
    end do
    call detail_leave(DP_OUTER_KERNEL)

  end subroutine

  subroutine begin_inverse203(oracle)
    logical, intent(in) :: oracle
    if(stage203_inverse.and.sparse_outer_install) call die("Stage203 requires sparse outer off")
    resident203=stage203_inverse
    oracle203=oracle
    work203=0_int64
  end subroutine

  subroutine end_inverse203
    resident203=.false.
  end subroutine

  subroutine check_gather203(plan,component,family)
    type(Transfer_Plan), intent(in) :: plan
    integer, intent(in) :: component,family
    integer :: lo,hi,p,k,v
    real(dp), allocatable :: expected(:,:,:)
    lo=merge(1,nscalar+1,component==1);hi=merge(nscalar,nscalar+EDGE,component==1)
    allocate(expected(lo:hi,inverse_levels,size(plan%slot)))
    do p=1,size(plan%slot)
       expected(:,:,p)=value(lo:hi,:,plan%slot(p),family)
    end do
    call transfer_nodes(plan,component,family,.false.)
    do p=1,size(plan%slot)
       do k=1,inverse_levels
          do v=lo,hi
             if (transfer_bits(expected(v,k,p))/=transfer_bits(value(v,k,plan%slot(p),family))) then
                write(error_unit,*) 'Stage203 phase mismatch rank/component/node/k/v ',rank,component,plan%slot(p),k,v
                write(error_unit,*) expected(v,k,p),value(v,k,plan%slot(p),family)
                call die('resident phase differs from independent block kernel')
             end if
          end do
       end do
    end do
  end subroutine

  subroutine remember_outer203(plan)
    type(Transfer_Plan), intent(in) :: plan
    integer :: p,q,s
    if (plan%edgewise) call die('Stage203 requires Stage192 sparse outer off')
    do p=1,size(plan%slot)
       s=plan%slot(p);q=restore203_index(s)
       if(q>0) restore203(nscalar+1:nscalar+EDGE,:,q)=value(nscalar+1:nscalar+EDGE,:,s,1)
    end do
  end subroutine

  subroutine store203(slot,component,k,x)
    integer, intent(in) :: slot,component,k
    real(dp), intent(in) :: x
    integer :: q
    value(component,k,slot,1)=x
    q=restore203_index(slot)
    if(q>0) restore203(component,k,q)=x
  end subroutine

  subroutine install_inverse203
    ! Install canonical interiors without exporting alias-overwritten values.
    ! Boundary views keep the latest dependency epoch; physics aliases have
    ! their independent consumer snapshot just as in Stage202.
    integer :: q,c
    if (.not.resident203) return
    do q=1,size(restore203_slot)
       saved203(:,:,q)=value(:,:,restore203_slot(q),1)
       value(:,:,restore203_slot(q),1)=restore203(:,:,q)
    end do
    if(ledger203) work203(6)=work203(6)+int(size(interior_plan%slot),int64)*(nscalar+EDGE)*inverse_levels
    do c=1,2
       call transfer_nodes(interior_plan,c,BLOCK_PAYLOAD_SOL,.true.)
    end do
    do q=1,size(restore203_slot)
       value(:,:,restore203_slot(q),1)=saved203(:,:,q)
    end do
    do c=1,2
       call transfer_nodes(boundary_plan,c,BLOCK_PAYLOAD_SOL,.true.)
    end do
  end subroutine

  subroutine compile203(total)
    integer, intent(in) :: total
    integer :: d,l,p,pi,c,ch,i,j,x,y,ip,jp,n,m,t,a,offp(N_BDRY+1),offc(N_BDRY+1)
    integer :: dimp(2,N_BDRY+1),dimc(2,N_BDRY+1),sx(11),sy(11),q,baseid
    integer :: tx(4),ty(4),wi(11)
    type(Scalar203) :: op
    type(Inner203) :: cell
    if (allocated(scalar203_ops)) deallocate(scalar203_ops,inner203_ops)
    allocate(scalar203_ops(max(1,2*total)),inner203_ops(max(1,total)))
    n=0;m=0
    wi=[1,2,3,1,2,3,4,1,2,3,4]
    do l=level_start,level_end-1
       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             if (.not.any(grid(d)%patch%elts(p+1)%children>0)) cycle
             call get_offs_Domain(grid(d),p,offp,dimp)
             do c=1,N_CHDRN
                ch=grid(d)%patch%elts(p+1)%children(c)
                if(ch==0) cycle
                call get_offs_Domain(grid(d),ch,offc,dimc)
                do j=1,PATCH_SIZE/2
                   y=2*(j-1);jp=j-1+chd_offs(2,c)
                   do i=1,PATCH_SIZE/2
                      x=2*(i-1);ip=i-1+chd_offs(1,c)
                      tx=[x,x+1,x+1,x];ty=[y,y,y+1,y+1]
                      do t=1,4
                         op=Scalar203(d=d,level=l,phase=merge(1,2,t==1))
                         op%target=node_slot(d,idx(tx(t),ty(t),offc,dimc))
                         op%mask=idx(x,y,offc,dimc)
                         if(t==1) then
                            op%n=11;sx=x+[1,1,2,0,-1,-1,-2,-1,0,-1,1]
                            sy=y+[0,1,1,1,0,1,-1,-1,-1,-2,-1]
                            baseid=idx(ip,jp,offp,dimp)
                            op%base=node_slot(d,baseid);op%area=baseid
                         else
                            op%n=4
                            select case(t)
                            case(2)
                               sx(:4)=[x,x+2,x+2,x];sy(:4)=[y,y,y+2,y-2]
                            case(3)
                               sx(:4)=[x+2,x,x+2,x];sy(:4)=[y+2,y,y,y+2]
                            case(4)
                               sx(:4)=[x,x,x-2,x+2];sy(:4)=[y,y+2,y,y+2]
                            end select
                            op%base=op%target;op%area=idx(tx(t),ty(t),offc,dimc)
                         end if
                         do q=1,op%n
                            a=idx(sx(q),sy(q),offc,dimc)
                            op%source(q)=node_slot(d,a)
                            op%geometry(q)=merge(a,op%area,t==1)
                         end do
                         n=n+1
                         if(n>size(scalar203_ops)) call die('scalar tape capacity')
                         scalar203_ops(n)=op
                      end do
                      cell=Inner203(d=d,level=l)
                      do t=LORT,UPLT
                         call edge_geom(4*(t-LORT)+1,ip,jp,DG,.true.)
                         call edge_geom(4*(t-LORT)+2,ip-t+1,jp,UP,.true.)
                         call edge_geom(4*(t-LORT)+3,ip,jp+t,RT,.true.)
                         call triangle_geom(4*(t-LORT)+4,ip,jp,t,.true.)
                      end do
                      call triangle_geom(9,x,y,LORT,.false.)
                      call edge_geom(10,x,y,RT,.false.)
                      call edge_geom(11,x,y,DG,.false.)
                      call edge_geom(12,x+1,y,UP,.false.)
                      call triangle_geom(13,x+1,y,LORT,.false.)
                      call edge_geom(14,x+1,y,RT,.false.)
                      call edge_geom(15,x+2,y,UP,.false.)
                      call edge_geom(16,x+1,y,DG,.false.)
                      call triangle_geom(17,x+1,y+1,LORT,.false.)
                      call edge_geom(18,x+1,y+1,DG,.false.)
                      call edge_geom(19,x+2,y+1,UP,.false.)
                      call edge_geom(20,x+1,y+1,RT,.false.)
                      call triangle_geom(21,x,y,UPLT,.false.)
                      call edge_geom(22,x,y,UP,.false.)
                      call edge_geom(23,x,y,DG,.false.)
                      call edge_geom(24,x,y+1,RT,.false.)
                      call triangle_geom(25,x,y+1,UPLT,.false.)
                      call edge_geom(26,x,y+1,UP,.false.)
                      call edge_geom(27,x,y+2,RT,.false.)
                      call edge_geom(28,x,y+1,DG,.false.)
                      call triangle_geom(29,x+1,y+1,UPLT,.false.)
                      call edge_geom(30,x+1,y+1,DG,.false.)
                      call edge_geom(31,x+1,y+2,RT,.false.)
                      call edge_geom(32,x+1,y+1,UP,.false.)
                      cell%target=cell%source([12,16,20,24,28,32])
                      m=m+1
                      if(m>size(inner203_ops)) call die('inner tape capacity')
                      inner203_ops(m)=cell
                   end do
                end do
             end do
          end do
       end do
    end do
    scalar203_ops=scalar203_ops(:n);inner203_ops=inner203_ops(:m)
  contains
    subroutine edge_geom(slot,xx,yy,e,parent)
      integer, intent(in) :: slot,xx,yy,e
      logical, intent(in) :: parent
      integer :: id
      if(parent) then
         id=idx(xx,yy,offp,dimp)
      else
         id=idx(xx,yy,offc,dimc)
      end if
      cell%node(slot)=id;cell%edge(slot)=e
      cell%source(slot)=edge_slot(d,EDGE*id+e)
    end subroutine
    subroutine triangle_geom(slot,xx,yy,e,parent)
      integer, intent(in) :: slot,xx,yy,e
      logical, intent(in) :: parent
      if(parent) then
         cell%node(slot)=idx(xx,yy,offp,dimp)
      else
         cell%node(slot)=idx(xx,yy,offc,dimc)
      end if
      cell%edge(slot)=e
    end subroutine
  end subroutine

  subroutine scalar_inverse203(level,phase)
    integer, intent(in) :: level,phase
    integer :: p,k,v,q,d
    integer, parameter :: wi(11)=[1,2,3,1,2,3,4,1,2,3,4]
    real(dp) :: w(11),area,s,x
    if(.not.resident203) return
    do p=1,size(scalar203_ops)
       associate(op=>scalar203_ops(p))
         if(op%level/=level.or.op%phase/=phase) cycle
         d=op%d
         if(phase==1) then
            if(grid(d)%mask_n%elts(op%mask+1)==FROZEN) cycle
         end if
         area=grid(d)%areas%elts(op%area+1)%hex_inv
         do q=1,op%n
            w(q)=grid(d)%overl_areas%elts(op%geometry(q)+1)%a(merge(wi(q),q,phase==1))
         end do
         do v=1,nscalar
            do k=1,inverse_levels
               s=0.0_dp
               do q=1,op%n
                  s=s+value(v,k,op%source(q),merge(2,1,phase==1))*w(q)
               end do
               if(phase==1) then
                  x=value(v,k,op%base,1)-s*area
               else
                  x=s*area+value(v,k,op%base,2)
               end if
               call store203(op%target,v,k,x)
            end do
         end do
         if(ledger203) work203(3)=work203(3)+int(inverse_levels,int64)*nscalar
       end associate
    end do
  end subroutine

  subroutine inner_inverse203(level)
    integer, intent(in) :: level
    integer :: p,k,q,t,d,id,e,s,slot,comp
    integer, parameter :: triangles(8)=[4,8,9,13,17,21,25,29]
    integer, parameter :: inputs(18)=[1,2,3,5,6,7,10,11,14,15,18,19,22,23,26,27,30,31]
    real(dp) :: g(32),u(32),curl(LORT:UPLT),x(6)
    if(.not.resident203) return
    do p=1,size(inner203_ops)
       associate(op=>inner203_ops(p))
         if(op%level/=level) cycle
         d=op%d
         do q=1,32
            id=op%node(q);e=op%edge(q)
            if(any(triangles==q)) then
               g(q)=grid(d)%triarea%elts(TRIAG*id+e+1)
            else
               g(q)=grid(d)%len%elts(EDGE*id+e+1)
            end if
         end do
         do k=1,inverse_levels
            do q=1,size(inputs)
               t=inputs(q);u(t)=alias_value(op%source(t),2,1,k,1)
            end do
            do t=LORT,UPLT
               s=4*(t-LORT)
               curl(t)=(u(s+1)*g(s+1)+u(s+2)*g(s+2)+u(s+3)*g(s+3))/g(s+4)
            end do
            x(1)=(g(9)*curl(LORT)-u(10)*g(10)-u(11)*g(11))/g(12)
            x(2)=(g(13)*curl(LORT)-u(14)*g(14)-u(15)*g(15))/g(16)
            x(3)=(g(17)*curl(LORT)-u(18)*g(18)-u(19)*g(19))/g(20)
            x(4)=(g(21)*curl(UPLT)-u(22)*g(22)-u(23)*g(23))/g(24)
            x(5)=(g(25)*curl(UPLT)-u(26)*g(26)-u(27)*g(27))/g(28)
            x(6)=(g(29)*curl(UPLT)-u(30)*g(30)-u(31)*g(31))/g(32)
            do q=1,6
               slot=(op%target(q)-1)/EDGE+1;comp=nscalar+mod(op%target(q)-1,EDGE)+1
               call store203(slot,comp,k,x(q)+value(comp,k,slot,2))
            end do
         end do
         if(ledger203) work203(4)=work203(4)+int(inverse_levels,int64)
       end associate
    end do
  end subroutine


  subroutine native_soil_transform
    ! Soil entry is the completed physics result, with its own retained aliases.
    ! Atmospheric coefficients are deliberately outside the selected transfer range.
    integer,allocatable::address(:,:)
    real(dp),allocatable::entry(:,:,:)
    integer(int64)::generation
    integer::l,d,k,v,q,c
    if(zmin>0) return
    call export_native_soil_aliases(address,entry,generation)
    call prepare_native_inverse(generation,first_level=level_start,profile=.false., &
         field_first=zmin,field_count=1-zmin,entry_alias_address=address,entry_alias_value=entry)
#ifdef SOIL_NATIVE_ORACLE
    call compare_native_inverse(sol(:,zmin:0),1,exact=.true.)
    call compare_native_inverse(sol(:,zmin:0),2,exact=.true.)
    write(*,'(a,i0)') 'Native soil entry exact rank=',rank
#endif
    do k=1,inverse_levels
       soil_layer=k
       do d=1,size(grid)
          soil_domain=d
          call apply_interscale_d(soil_Restrict_velo,grid(d),level_start-1,k,0,0)
       end do
    end do
    do q=1,size(scaffold_slot)
       scaffold_value(:,:,q)=value(:,:,scaffold_slot(q),1)
    end do
    do c=1,2
       call native_inverse_boundary(c,BLOCK_PAYLOAD_SOL,level_start-1,level_end)
    end do
    do l=level_start-1,level_end-1
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             do v=1,nscalar
                soil_variable=v
                call apply_interscale_d(soil_Compute_scalar_wavelets,grid(d),l,Z_NULL,0,0)
             end do
             call apply_interscale_d(soil_Compute_velo_wavelets,grid(d),l,Z_NULL,0,0)
          end do
       end do
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             call apply_to_penta_d(soil_Compute_velo_wavelets_penta,grid(d),l,Z_NULL)
          end do
       end do
    end do
    do c=1,2
       call native_inverse_boundary(c,BLOCK_PAYLOAD_WAV_COEFF,level_start-1,level_end)
    end do
    do k=1,inverse_levels
       soil_layer=k
       do l=level_start+1,level_end
          do d=1,size(grid)
             soil_domain=d
             do v=1,nscalar
                soil_variable=v
                call apply_onescale_d(soil_compress_scalar,grid(d),l,Z_NULL,0,1)
             end do
             call apply_onescale_d(soil_compress_vector,grid(d),l,Z_NULL,0,0)
          end do
       end do
    end do
    do c=1,2
       call transfer_nodes(interior_plan,c,BLOCK_PAYLOAD_WAV_COEFF,.true.)
       call native_inverse_boundary(c,BLOCK_PAYLOAD_WAV_COEFF,level_start,level_end)
       call native_inverse_boundary(c,BLOCK_PAYLOAD_SOL,level_start,level_end)
    end do
    ! Unlike the canonical atmospheric tape, soil reconstruction also writes
    ! geometric aliases. Keep the complete legacy walk on native operands.
    do l=level_start,level_end-1
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             do v=1,nscalar
                soil_variable=v
                call apply_interscale_d2(soil_Prolong_scalar,grid(d),l,Z_NULL,0,1)
             end do
          end do
       end do
       call native_inverse_boundary(1,BLOCK_PAYLOAD_SOL,l+1,l+1)
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             call apply_interscale_d2(soil_Reconstruct_outer_velo,grid(d),l,Z_NULL,0,1)
             call apply_to_penta_d(soil_Reconstruct_velo_penta,grid(d),l,Z_NULL)
          end do
       end do
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             do v=1,nscalar
                soil_variable=v
                call apply_interscale_d(soil_Reconstruct_scalar,grid(d),l,Z_NULL,0,0)
             end do
          end do
       end do
       call native_inverse_boundary(2,BLOCK_PAYLOAD_SOL,l+1,l+1)
       do k=1,inverse_levels
          soil_layer=k
          do d=1,size(grid)
             soil_domain=d
             call apply_interscale_d(soil_Reconstruct_inner_velo,grid(d),l,Z_NULL,0,0)
          end do
       end do
       if(l<level_end-1) call native_inverse_boundary(2,BLOCK_PAYLOAD_SOL,l+1,l+1)
    end do
    do q=1,size(consumer_alias_slot)
       consumer_alias_value(:,:,q)=value(:,:,consumer_alias_slot(q),1)
    end do
    do c=1,2
       call native_inverse_scatter(c,BLOCK_PAYLOAD_SOL,.true.)
    end do
    write(*,'(a,i0,a,i0)') 'Native soil transform rank=',rank,' layers=',inverse_levels
  end subroutine native_soil_transform



  subroutine publish_native_soil(scaling,wavelet)
    type(Float_Field),intent(inout)::scaling(:,:),wavelet(:,:)
    if(inverse_field_first/=zmin.or.inverse_levels/=1-zmin) call die('soil publication epoch')
    call publish_native_inverse(scaling,all_nodes=.true.)
    call publish_native_inverse(scaling,consumer_aliases=.true.)
    call publish_native_inverse(wavelet,payload_family=BLOCK_PAYLOAD_WAV_COEFF,all_nodes=.true.)
    scaling%bdry_uptodate=.false.
    wavelet%bdry_uptodate=.false.
  end subroutine publish_native_soil

  function soil_scalar(address) result(x)
    integer,intent(in)::address
    real(dp),pointer::x
    integer::id,p,v
    id=address-1
    v=soil_variable
    if(soil_domain<1.or.soil_domain>size(nodes)) call die('soil callback Domain')
    if(id<0.or.id>=size(nodes(soil_domain)%slot)) call die('soil callback address')
    p=nodes(soil_domain)%slot(id+1)
    if(p<1) then
       write(error_unit,*) 'Missing soil operand rank/domain/id ',rank,soil_domain,id
       call die('soil dependency not registered')
    end if
    x=>value(v,soil_layer,p,1)
  end function soil_scalar

  function soil_vector(address) result(x)
    integer,intent(in)::address
    real(dp),pointer::x
    integer::id,p,v
    id=(address-1)/EDGE
    v=nscalar+mod(address-1,EDGE)+1
    if(soil_domain<1.or.soil_domain>size(nodes)) call die('soil callback Domain')
    if(id<0.or.id>=size(nodes(soil_domain)%slot)) call die('soil callback address')
    p=nodes(soil_domain)%slot(id+1)
    if(p<1) then
       write(error_unit,*) 'Missing soil operand rank/domain/id ',rank,soil_domain,id
       call die('soil dependency not registered')
    end if
    x=>value(v,soil_layer,p,1)
  end function soil_vector

  function soil_scalar_wavelet(address) result(x)
    integer,intent(in)::address
    real(dp),pointer::x
    integer::id,p,v
    id=address-1
    v=soil_variable
    if(soil_domain<1.or.soil_domain>size(nodes)) call die('soil callback Domain')
    if(id<0.or.id>=size(nodes(soil_domain)%slot)) call die('soil callback address')
    p=nodes(soil_domain)%slot(id+1)
    if(p<1) then
       write(error_unit,*) 'Missing soil operand rank/domain/id ',rank,soil_domain,id
       call die('soil dependency not registered')
    end if
    x=>value(v,soil_layer,p,2)
  end function soil_scalar_wavelet

  function soil_vector_wavelet(address) result(x)
    integer,intent(in)::address
    real(dp),pointer::x
    integer::id,p,v
    id=(address-1)/EDGE
    v=nscalar+mod(address-1,EDGE)+1
    if(soil_domain<1.or.soil_domain>size(nodes)) call die('soil callback Domain')
    if(id<0.or.id>=size(nodes(soil_domain)%slot)) call die('soil callback address')
    p=nodes(soil_domain)%slot(id+1)
    if(p<1) then
       write(error_unit,*) 'Missing soil operand rank/domain/id ',rank,soil_domain,id
       call die('soil dependency not registered')
    end if
    x=>value(v,soil_layer,p,2)
  end function soil_vector_wavelet

  subroutine set_soil_scalar_wavelet(address,x)
    integer,intent(in)::address
    real(dp),intent(in)::x
    real(dp),pointer::destination
    destination=>soil_scalar_wavelet(address)
    destination=x
  end subroutine set_soil_scalar_wavelet

  subroutine set_soil_vector_wavelet(address,x)
    integer,intent(in)::address
    real(dp),intent(in)::x
    real(dp),pointer::destination
    destination=>soil_vector_wavelet(address)
    destination=x
  end subroutine set_soil_vector_wavelet

  subroutine set_soil_scalar(address,x)
    integer,intent(in)::address
    real(dp),intent(in)::x
    real(dp),pointer::destination
    destination=>soil_scalar(address)
    destination=x
  end subroutine set_soil_scalar

  subroutine set_soil_vector(address,x)
    integer,intent(in)::address
    real(dp),intent(in)::x
    real(dp),pointer::destination
    destination=>soil_vector(address)
    destination=x
  end subroutine set_soil_vector

  subroutine soil_Compute_scalar_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute wavelet coefficients for scalars
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd    = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)

    if (dom%mask_n%elts(idE_chd+1) >= ADJZONE)  call set_soil_scalar_wavelet(idE_chd+1,soil_scalar(idE_chd+1) - &
         & soil_Interp_node (dom, idE_chd, id_chd, id2E_chd,  id2NE_chd, id2S_chd))

    if (dom%mask_n%elts(idNE_chd+1) >= ADJZONE)  call set_soil_scalar_wavelet(idNE_chd+1,soil_scalar(idNE_chd+1) - &
         & soil_Interp_node (dom, idNE_chd, id2NE_chd,  id_chd, id2E_chd, id2N_chd))

    if (dom%mask_n%elts(idN_chd+1) >= ADJZONE)  call set_soil_scalar_wavelet(idN_chd+1,soil_scalar(idN_chd+1) - &
         & soil_Interp_node (dom, idN_chd, id_chd, id2N_chd,  id2W_chd, id2NE_chd))
  end subroutine soil_Compute_scalar_wavelets
  subroutine soil_Compute_velo_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute velocity wavelet coefficients (except at pentagons)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer                :: e, idN_chd, idE_chd, idNE_chd, id1, id2
    real(dp)               :: u
    real(dp), dimension(6) :: u_inner

    ! Velocity wavelets at 6 outer child edges
    do e = 1, EDGE
       id1 = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2 = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       if (dom%mask_e%elts(EDGE*id2+e) < ADJZONE) cycle

       u = soil_Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)

       call set_soil_vector_wavelet(EDGE*id2+e,soil_vector(EDGE*id2+e) - u)
       call set_soil_vector_wavelet(EDGE*id1+e,- soil_vector_wavelet(EDGE*id2+e)) ! Restriction(Prolongation) = Identity
    end do

    ! Velocity wavelets at 6 inner child edges
    u_inner = soil_Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*idE_chd+UP+1)  >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idE_chd+UP+1,soil_vector(EDGE*idE_chd+UP+1)  - u_inner(1))
    if (dom%mask_e%elts(EDGE*idE_chd+DG+1)  >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idE_chd+DG+1,soil_vector(EDGE*idE_chd+DG+1)  - u_inner(2))
    if (dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idNE_chd+RT+1,soil_vector(EDGE*idNE_chd+RT+1) - u_inner(3))
    if (dom%mask_e%elts(EDGE*idN_chd+RT+1)  >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idN_chd+RT+1,soil_vector(EDGE*idN_chd+RT+1)  - u_inner(4))
    if (dom%mask_e%elts(EDGE*idN_chd+DG+1)  >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idN_chd+DG+1,soil_vector(EDGE*idN_chd+DG+1)  - u_inner(5))
    if (dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE) call &
         & set_soil_vector_wavelet(EDGE*idNE_chd+UP+1,soil_vector(EDGE*idNE_chd+UP+1) - u_inner(6))
  end subroutine soil_Compute_velo_wavelets
  subroutine soil_Compute_velo_wavelets_penta (dom, p, c, offs, dims, zlev)
    ! Compute velocity wavelet coefficients at pentagons.
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, c, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id_chd, idE_chd, idN_chd, p_chd
    integer :: offs_chd(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1)

    real(dp) :: v
    real(dp) :: u(2)

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then
       ! Parts 3 and 4 of hexagon IMINUSJPLUS, at the upper-left
       ! corner of the lozenge, combine to form a pentagon.
       !
       ! pedlen(EDGE*idW + RT + 1) = 0 in this case.
       id_chd  = idx(0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx(0, LAST,   offs_chd, dims_chd)

       v = -(  base_weight(8) +  real(dom%I_u_wgt%elts(idN_chd+1)%enc(8), kind=dp)  ) * (  soil_vector(idx(0,  &
         & PATCH_SIZE, offs, dims)*EDGE + UP + 1) +  soil_vector(idx(-1, PATCH_SIZE, offs, dims)*EDGE + RT + 1)    )

       if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
          call set_soil_vector_wavelet(EDGE*id_chd +UP+1,soil_vector_wavelet(EDGE*id_chd +UP+1) - v)

          call set_soil_vector_wavelet(EDGE*idN_chd+UP+1,soil_vector_wavelet(EDGE*idN_chd+UP+1) + v)
       end if

    else if (c == IPLUSJMINUS) then
       ! Parts 5 and 6 of hexagon IPLUSJMINUS, at the lower-right
       ! corner of the lozenge, combine to form a pentagon.
       !
       ! pedlen(EDGE*idS + UP + 1) = 0 in this case.
       id_chd  = idx(LAST-1, 0, offs_chd, dims_chd)
       idE_chd = idx(LAST,   0, offs_chd, dims_chd)

       v = (  base_weight(7) +  real(dom%I_u_wgt%elts(idE_chd+1)%enc(7), kind=dp)  ) * (  soil_vector(idx(PATCH_SIZE,  &
         & 0, offs, dims)*EDGE + RT + 1) +  soil_vector(idx(PATCH_SIZE, -1, offs, dims)*EDGE + UP + 1)    )

       if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
          call set_soil_vector_wavelet(EDGE*id_chd +RT+1,soil_vector_wavelet(EDGE*id_chd +RT+1) - v)

          call set_soil_vector_wavelet(EDGE*idE_chd+RT+1,soil_vector_wavelet(EDGE*idE_chd+RT+1) + v)
       end if
    end if

    if (c /= IJMINUS) return

    ! Parts 4 and 5 of hexagon IJMINUS, at the lower-left corner
    ! of the lozenge, combine to form a pentagon.
    !
    ! pedlen(EDGE*idSW + DG + 1) = 0 in this case.
    id_chd  = idx(0, 0, offs_chd, dims_chd)
    idN_chd = idx(0, 1, offs_chd, dims_chd)
    idE_chd = idx(1, 0, offs_chd, dims_chd)

    u = soil_velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
       call set_soil_vector_wavelet(EDGE*id_chd +UP+1,soil_vector_wavelet(EDGE*id_chd +UP+1) + u(1))

       call set_soil_vector_wavelet(EDGE*idN_chd+UP+1,soil_vector_wavelet(EDGE*idN_chd+UP+1) - u(1))
    end if

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
       call set_soil_vector_wavelet(EDGE*id_chd +RT+1,soil_vector_wavelet(EDGE*id_chd +RT+1) + u(2))

       call set_soil_vector_wavelet(EDGE*idE_chd+RT+1,soil_vector_wavelet(EDGE*idE_chd+RT+1) - u(2))
    end if
  end subroutine soil_Compute_velo_wavelets_penta
  function soil_velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd) result(val)

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer, dimension(N_BDRY+1),   intent(in)    :: offs, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims, dims_chd
    real(dp)                                      :: val(2)

    integer :: i, j, i_chd, j_chd

    i     = 0
    j     = 0
    i_chd = 0
    j_chd = 0

    val = [   (base_weight(9)  + real (dom%I_u_wgt%elts(idx__fast(i_chd+end_pt(1,2,UP+1), j_chd+end_pt(2,2,UP+1), &
         & offs_chd(1))+1)%enc(9),kind=dp)) *  ((-soil_vector(idx(0, -1, offs, dims)*EDGE+UP+1) - &
         & (-soil_vector(idx(-1, -1, offs, dims)*EDGE+1))) -  (soil_vector(ed_idx(i+end_pt(1,1,UP+1), &
         & j+end_pt(2,1,UP+1), hex_sides(:,hex_s_offs(UP+1)+0+1), offs, dims)+1) -  &
         & soil_vector(ed_idx(i+opp_no(1,2,UP+1), j+opp_no(2,2,UP+1), hex_sides(:,hex_s_offs(UP+1)+1+1), offs, &
         & dims)+1))),  (base_weight(6)  + real (dom%I_u_wgt%elts(idx__fast(i_chd + end_pt(1,2,RT+1), j_chd + &
         & end_pt(2,2,RT+1), offs_chd(1))+1)%enc(6),kind=dp)) *  (soil_vector(idx(-1, -1, offs, dims)*EDGE+1) + &
         & soil_vector(idx(-1, 0, offs, dims)*EDGE + RT+1)  - (soil_vector(ed_idx(i+opp_no(1,1,RT+1), &
         & j+opp_no(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+1+1), offs, dims)+1)  -  &
         & soil_vector(ed_idx(i+end_pt(1,1,RT+1), j+end_pt(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+2+1), offs, dims)+1))) ]
  end function soil_velo_interp_penta_corr
  function soil_Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd) result(val)
    ! Interpolate inner velocities to fine edges

    implicit none

    type(Domain), intent (in) :: dom
    integer,      intent(in)  :: i_par, j_par, i_chd, j_chd
    integer,      intent(in)  :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)  :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)
    real(dp)                  :: val(6)

    integer                        :: id, id_par, id1_par, id2_par, t, idN, idE, idNE, idN2E, id2NE, idN2, idE2
    real(dp), dimension(LORT:UPLT) :: curl_u

    id_par = idx(i_par, j_par, offs_par, dims_par)

    curl_u = 0.0_dp
    do t = LORT, UPLT
       id1_par = idx (i_par-t+1, j_par,   offs_par, dims_par)
       id2_par = idx (i_par,     j_par+t, offs_par, dims_par)
       curl_u(t) = (soil_vector(EDGE*id_par +DG+1)*dom%len%elts(EDGE*id_par+DG+1)  + &
         & soil_vector(EDGE*id1_par+UP+1)*dom%len%elts(EDGE*id1_par+UP+1)  + &
         & soil_vector(EDGE*id2_par+RT+1)*dom%len%elts(EDGE*id2_par+RT+1)) / dom%triarea%elts(TRIAG*id_par+t+1)
    end do

    id = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN2  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    idE2  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)

    val(1) = (dom%triarea%elts(TRIAG*id+LORT+1) * curl_u(LORT)  - soil_vector(EDGE*id+RT+1) * &
         & dom%len%elts(EDGE*id+RT+1) - soil_vector(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1))  / dom%len%elts(EDGE*idE+UP+1)

    val(2) = (dom%triarea%elts(TRIAG*idE+LORT+1) * curl_u(LORT)  - soil_vector(EDGE*idE+RT+1) * &
         & dom%len%elts(EDGE*idE+RT+1) - soil_vector(EDGE*idE2+UP+1) * dom%len%elts(EDGE*idE2+UP+1))  / dom%len%elts(EDGE*idE+DG+1)

    val(3) = (dom%triarea%elts(TRIAG*idNE+LORT+1) * curl_u(LORT)  - soil_vector(EDGE*idNE+DG+1) * &
         & dom%len%elts(EDGE*idNE+DG+1) - soil_vector(EDGE*idN2E+UP+1) *  dom%len%elts(EDGE*idN2E+UP+1))  / &
         & dom%len%elts(EDGE*idNE+RT+1)

    val(4) = (dom%triarea%elts(TRIAG*id+UPLT+1) * curl_u(UPLT)  - soil_vector(EDGE*id+UP+1) * &
         & dom%len%elts(EDGE*id+UP+1) - soil_vector(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1))  / dom%len%elts(EDGE*idN+RT+1)

    val(5) = (dom%triarea%elts(TRIAG*idN+UPLT+1) * curl_u(UPLT)  - soil_vector(EDGE*idN+UP+1) * &
         & dom%len%elts(EDGE*idN+UP+1) - soil_vector(EDGE*idN2+RT+1) * dom%len%elts(EDGE*idN2+RT+1))  / dom%len%elts(EDGE*idN+DG+1)

    val(6) = (dom%triarea%elts(TRIAG*idNE+UPLT+1) * curl_u(UPLT)  - soil_vector(EDGE*idNE+DG+1) * &
         & dom%len%elts(EDGE*idNE+DG+1) - soil_vector(EDGE*id2NE+RT+1) *  dom%len%elts(EDGE*id2NE+RT+1))  / &
         & dom%len%elts(EDGE*idNE+UP+1)
  end function soil_Interp_inner_velo
  function soil_Interp_outer_velo (dom, i_par, j_par, ide, e, offs, dims) result(val)
    ! Interpolate outer velocities to fine edge EDGE*ide + e
    implicit none

    type(Domain),                   intent(in) :: dom
    integer,                        intent(in) :: i_par, j_par, ide, e
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims
    real(dp)                                   :: val

    real(dp), dimension(9) :: wgt

    wgt = base_weight + real (dom%I_u_wgt%elts(ide+1)%enc, kind=dp)

   val = sum (wgt * [ &
         soil_vector(EDGE * idx(i_par, j_par, offs, dims) + e), &
         soil_vector(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+2+1), offs, dims) + 1), &
         soil_vector(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+3+1), offs, dims) + 1), &
         soil_vector(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+5+1), offs, dims) + 1), &
         soil_vector(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+0+1), offs, dims) + 1), &

         soil_vector(ed_idx(i_par + opp_no(1,1,e), j_par + opp_no(2,1,e), hex_sides(:,hex_s_offs(e)+1+1), offs, dims) + 1) - &
         soil_vector(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+2+1), offs, dims) + 1), &

         soil_vector(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+3+1), offs, dims) + 1) - &
         soil_vector(ed_idx(i_par + opp_no(1,1,e), j_par + opp_no(2,1,e), hex_sides(:,hex_s_offs(e)+4+1), offs, dims) + 1), &

         soil_vector(ed_idx(i_par + opp_no(1,2,e), j_par + opp_no(2,2,e), hex_sides(:,hex_s_offs(e)+4+1), offs, dims) + 1) - &
         soil_vector(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+5+1), offs, dims) + 1), &

         soil_vector(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+0+1), offs, dims) + 1) - &
         soil_vector(ed_idx(i_par + opp_no(1,2,e), j_par + opp_no(2,2,e), hex_sides(:,hex_s_offs(e)+1+1), offs, dims) + 1) ])
  end function soil_Interp_outer_velo
  subroutine soil_Restrict_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Circulation conserving velocity  restriction (note that fine edges are exact bisections of coarse edges)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id_par

    id_par   = idx (i_par, j_par, offs_par, dims_par)
    id_chd   = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) > 0) call set_soil_vector(EDGE*id_par+RT+1,0.5 * &
         & (soil_vector(EDGE*id_chd+RT+1)    + soil_vector(EDGE*idE_chd+RT+1)))
    if (dom%mask_e%elts(EDGE*id_chd+DG+1) > 0) call set_soil_vector(EDGE*id_par+DG+1,0.5 * &
         & (soil_vector(EDGE*idNE_chd+DG+1)  + soil_vector(EDGE*id_chd+DG+1)))
    if (dom%mask_e%elts(EDGE*id_chd+UP+1) > 0) call set_soil_vector(EDGE*id_par+UP+1,0.5 * &
         & (soil_vector(EDGE*id_chd+UP+1)    + soil_vector(EDGE*idN_chd+UP+1)))
  end subroutine soil_Restrict_velo
  real(dp) function soil_Interp_node (dom, id, id1, id2, id3, id4)
    ! Interpolation at nodes
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: id, id1, id2, id3, id4

    soil_Interp_node = (  dom%overl_areas%elts(id+1)%a(1) * soil_scalar(id1+1) +  dom%overl_areas%elts(id+1)%a(2) * &
         & soil_scalar(id2+1) +  dom%overl_areas%elts(id+1)%a(3) * soil_scalar(id3+1) +  &
         & dom%overl_areas%elts(id+1)%a(4) * soil_scalar(id4+1)  ) * dom%areas%elts(id+1)%hex_inv
  end function soil_Interp_node
  subroutine soil_compress_scalar (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id_i) < ADJZONE) call set_soil_scalar_wavelet(id_i,0.0_dp)
  end subroutine soil_compress_scalar
  subroutine soil_compress_vector (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: e, id

    id = idx (i, j, offs, dims)

    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) < ADJZONE) call set_soil_vector_wavelet(EDGE*id+e,0.0_dp)
    end do
  end subroutine soil_compress_vector

  subroutine soil_Prolong_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
    integer :: id_par, id_chd, idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE
    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    if (dom%mask_n%elts(id_chd+1) == FROZEN) return
    idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idS2W = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)
    idSW  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    id2SW = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idSE  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)
    call set_soil_scalar(id_chd+1,soil_scalar(id_par+1) - (  soil_scalar_wavelet(idE+1)   * dom%overl_areas%elts(idE  &
         & +1)%a(1) +  soil_scalar_wavelet(idNE+1)  * dom%overl_areas%elts(idNE +1)%a(2) +  &
         & soil_scalar_wavelet(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) +  soil_scalar_wavelet(idN+1)   * &
         & dom%overl_areas%elts(idN  +1)%a(1) +  soil_scalar_wavelet(idW+1)   * dom%overl_areas%elts(idW  +1)%a(2) +  &
         & soil_scalar_wavelet(idNW+1)  * dom%overl_areas%elts(idNW +1)%a(3) +  soil_scalar_wavelet(idS2W+1) * &
         & dom%overl_areas%elts(idS2W+1)%a(4) +  soil_scalar_wavelet(idSW+1)  * dom%overl_areas%elts(idSW+ 1)%a(1) +  &
         & soil_scalar_wavelet(idS+1)   * dom%overl_areas%elts(idS  +1)%a(2) +  soil_scalar_wavelet(id2SW+1) * &
         & dom%overl_areas%elts(id2SW+1)%a(3) +  soil_scalar_wavelet(idSE+1)  * dom%overl_areas%elts(idSE +1)%a(4)  ) &
         & * dom%areas%elts(id_par+1)%hex_inv)
  end subroutine soil_Prolong_scalar

  subroutine soil_Reconstruct_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)
    call set_soil_scalar(idE_chd+1,soil_Interp_node (dom, idE_chd,  id_chd,    id2E_chd, id2NE_chd, id2S_chd)  + &
         & soil_scalar_wavelet(idE_chd+1))
    call set_soil_scalar(idNE_chd+1,soil_Interp_node (dom, idNE_chd, id2NE_chd, id_chd,   id2E_chd,  id2N_chd)  + &
         & soil_scalar_wavelet(idNE_chd+1))
    call set_soil_scalar(idN_chd+1,soil_Interp_node (dom, idN_chd,  id_chd,    id2N_chd, id2W_chd,  id2NE_chd) + &
         & soil_scalar_wavelet(idN_chd+1))
  end subroutine soil_Reconstruct_scalar

  subroutine soil_Reconstruct_outer_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
    integer  :: e, id_par, id1, id2
    real(dp) :: u
    id_par = idx (i_par, j_par, offs_par, dims_par)
    do e = 1, EDGE
       id1    = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2    = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)
       u = soil_Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)
       call set_soil_vector(EDGE*id2+e,u + soil_vector_wavelet(EDGE*id2+e))
       call set_soil_vector(EDGE*id1+e,2 * soil_vector(EDGE*id_par+e) - soil_vector(EDGE*id2+e))
    end do
  end subroutine soil_Reconstruct_outer_velo

  subroutine soil_Reconstruct_inner_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
    integer                :: id_chd, idN_chd, idE_chd, idNE_chd
    real(dp), dimension(6) :: u_inner
    id_chd   = idx (i_chd, j_chd, offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    u_inner = soil_Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)
    call set_soil_vector(EDGE*idE_chd +UP+1,u_inner(1) + soil_vector_wavelet(EDGE*idE_chd +UP+1))
    call set_soil_vector(EDGE*idE_chd +DG+1,u_inner(2) + soil_vector_wavelet(EDGE*idE_chd +DG+1))
    call set_soil_vector(EDGE*idNE_chd+RT+1,u_inner(3) + soil_vector_wavelet(EDGE*idNE_chd+RT+1))
    call set_soil_vector(EDGE*idN_chd +RT+1,u_inner(4) + soil_vector_wavelet(EDGE*idN_chd +RT+1))
    call set_soil_vector(EDGE*idN_chd +DG+1,u_inner(5) + soil_vector_wavelet(EDGE*idN_chd +DG+1))
    call set_soil_vector(EDGE*idNE_chd+UP+1,u_inner(6) + soil_vector_wavelet(EDGE*idNE_chd+UP+1))
  end subroutine soil_Reconstruct_inner_velo

  subroutine soil_Reconstruct_velo_penta (dom, p, c, offs, dims, z_null)
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: p, c, z_null
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    integer                        :: id_chd,  idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(dp)                       :: v
    real(dp), dimension(2)         :: u
    p_chd = dom%patch%elts(p+1)%children(c-4)
    if (p_chd == 0) return
    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)
    if (c == IMINUSJPLUS) then
       id_chd  = idx (0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx (0, LAST,   offs_chd, dims_chd)
       v = (base_weight(8) + real (dom%I_u_wgt%elts(idN_chd+1)%enc(8),kind=dp))*(  soil_vector(EDGE*idx( 0, &
         & PATCH_SIZE, offs, dims)+UP+1)  + soil_vector(EDGE*idx(-1, PATCH_SIZE, offs, dims)+RT+1))
       call set_soil_vector(EDGE*id_chd+UP+1,soil_vector(EDGE*id_chd +UP+1) - v)
       call set_soil_vector(EDGE*idN_chd+UP+1,soil_vector(EDGE*idN_chd+UP+1) + v)
    else
       if (c == IPLUSJMINUS) then
          id_chd  = idx (LAST-1, 0, offs_chd, dims_chd)
          idE_chd = idx (LAST,   0, offs_chd, dims_chd)
          v = -(base_weight(7) + real (dom%I_u_wgt%elts(idE_chd+1)%enc(7),kind=dp))*(  &
         & soil_vector(EDGE*idx(PATCH_SIZE,  0, offs, dims)+RT+1)  + soil_vector(EDGE*idx(PATCH_SIZE, -1, offs, dims)+UP+1))
          call set_soil_vector(EDGE*id_chd +RT+1,soil_vector(EDGE*id_chd +RT+1) - v)
          call set_soil_vector(EDGE*idE_chd+RT+1,soil_vector(EDGE*idE_chd+RT+1) + v)
       end if
    end if
    if (.not. c == IJMINUS) return
    id_chd  = idx (0, 0, offs_chd, dims_chd)
    idN_chd = idx (0, 1, offs_chd, dims_chd)
    idE_chd = idx (1, 0, offs_chd, dims_chd)
    u = soil_velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)
    call set_soil_vector(EDGE*id_chd +UP+1,soil_vector(EDGE*id_chd +UP+1) - u(1))
    call set_soil_vector(EDGE*idN_chd+UP+1,soil_vector(EDGE*idN_chd+UP+1) + u(1))
    call set_soil_vector(EDGE*id_chd +RT+1,soil_vector(EDGE*id_chd +RT+1) - u(2))
    call set_soil_vector(EDGE*idE_chd+RT+1,soil_vector(EDGE*idE_chd+RT+1) + u(2))
  end subroutine soil_Reconstruct_velo_penta
end module parallel_block_inverse_mod
