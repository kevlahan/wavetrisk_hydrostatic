module parallel_block_columns_mod
  ! Interstage hydrostatic producer. Compact interiors are read directly on
  ! their final block owners. Only explicitly requested inverse aliases are
  ! delivered as operands. Geometry owners receive results for the remaining
  ! primitive/physics consumers, never a replacement solution field.
  use iso_fortran_env, only: int64
  use ieee_arithmetic, only: ieee_is_finite
  use kind_mod, only: dp
  use mpi_f08
  use shared_mod, only: N_BDRY,S_MASS,S_TEMP,zlevels,level_start,level_end,N_CHDRN, &
       grav_accel,p_top,c_p,p_0,kappa,compressible
  use arch_mod, only: comm,rank,n_process,block_catalog,owner,glo_id,abort_run
  use domain_mod, only: Domain,grid,idx,sol_mean
  use domain_ops_mod, only: apply_onescale_to_patch
  use parallel_block_mod, only: Block_Data,STORE_PATCH,apply_local_block_field_consumer
  use parallel_block_inverse_mod, only: inverse_operand_key,inverse_scalar_operand,inverse_operand_generation
  use parallel_block_inverse_mod, only: require_tendency_operands
  use parallel_block_primitive_kernels_mod, only: pressure_layer
  use init_mod, only: surf_geopot
  implicit none
  private
  public :: arm_block_columns,start_block_columns,column_surface,begin_column_phase,column_result,end_block_columns
  public :: columns_active,columns_oracle,columns_without_reference
  public :: owner_tendency_without_reference
  public :: columns_initial_pending
  logical :: initial_pending=.false.
  public :: column_consumer_address,owner_column_result
  logical :: armed=.false.,columns_active=.false.,columns_oracle=.false.
  type :: Node_Index
     integer, allocatable :: column(:)
  end type
  type :: Column_Context
     integer :: layer=0,level=0
  end type
  type(Node_Index), allocatable :: lookup(:)
  integer, allocatable :: local_key(:,:),local_domain(:),local_node(:),visits(:,:),surface_visits(:)
  integer, allocatable :: sc(:),sd(:),rc(:),rd(:),request(:,:),order(:),packed_local(:),block_head(:),next_column(:)
  integer, allocatable :: alias_sc(:),alias_sd(:),alias_rc(:),alias_rd(:),alias_index(:)
  integer, allocatable :: phase_sc(:),phase_sd(:),phase_rc(:),phase_rd(:),first(:),cursor(:),remote_first(:)
  real(dp), allocatable :: alias_operand(:,:,:),lower(:),geopot(:),surface(:),send_result(:,:),result(:,:)
  logical, allocatable :: filled(:)
  integer(int64) :: generation=-1_int64,executed=0_int64,checked=0_int64
  integer :: ncolumn=0,compile_level=0,phase_layer=0,phase_level=0,expected_layer=1,expected_level=0
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': block columns: ',message
    call abort_run
    error stop 'block columns abort returned'
  end subroutine

  subroutine check(ierr,where)
    integer, intent(in) :: ierr
    character(*), intent(in) :: where
    if(ierr/=MPI_SUCCESS) call die(where)
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

  subroutine arm_block_columns(stage,oracle)
    integer, intent(in) :: stage
    logical, intent(in) :: oracle
    if(columns_active.or.armed) call die('unfinished producer transaction')
#ifdef BLOCK_COLUMNS
    ! Stage one can use the distinct final-transform alias epoch and remains
    ! on the compatibility producer until that lifecycle is migrated.
    armed=stage>1
#endif
#ifdef BLOCK_TENDENCY_INITIAL
    armed=stage>=1;initial_pending=stage==1
#endif
    columns_oracle=oracle
#ifdef BLOCK_COLUMNS_ORACLE
    columns_oracle=.true.
#endif
  end subroutine

  logical function columns_initial_pending()
    columns_initial_pending=armed.and.initial_pending
  end function

  logical function owner_tendency_without_reference()
    ! One predicate for removing compatibility products. Any independent
    ! reference selector keeps those products available to its old kernels.
    owner_tendency_without_reference=.false.
#if defined(BLOCK_VELOCITY_OWNER) && !defined(BLOCK_COLUMNS_ORACLE) && !defined(BLOCK_PRIMITIVES_ORACLE) && !defined(BLOCK_DIFFUSION_ORACLE) && !defined(BLOCK_MASS_OWNER_ORACLE) && !defined(BLOCK_VELOCITY_OWNER_ORACLE)
    owner_tendency_without_reference=columns_active.and..not.columns_oracle
#endif
  end function

  logical function columns_without_reference()
    columns_without_reference=armed.and..not.columns_oracle
  end function

  subroutine register_visit(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: d,id,c
    d=dom%id+1;id=idx(i,j,offs,dims)+1
    c=lookup(d)%column(id)
    if(c==0) then
       ncolumn=ncolumn+1;c=ncolumn
       lookup(d)%column(id)=c
       local_domain(c)=d;local_node(c)=id
    end if
    if(compile_level==0) then
       surface_visits(c)=surface_visits(c)+1
    else
       visits(compile_level,c)=visits(compile_level,c)+1
    end if
  end subroutine

  subroutine compile_columns
    integer :: n,d,p,l,j,c,r,b,pos,ierr,h,a
    integer, allocatable :: send_key(:,:),counts(:)
    if(allocated(lookup)) then
       deallocate(lookup,local_key,local_domain,local_node,visits,surface_visits,sc,sd,rc,rd, &
            request,order,packed_local,block_head,next_column,alias_sc,alias_sd,alias_rc,alias_rd,alias_index, &
            phase_sc,phase_sd,phase_rc,phase_rd,first,cursor,remote_first,lower,geopot,surface)
    end if
    allocate(lookup(size(grid)))
    n=0
    do d=1,size(grid)
       allocate(lookup(d)%column(grid(d)%node%length));lookup(d)%column=0
       n=n+grid(d)%node%length
    end do
    allocate(local_domain(n),local_node(n),visits(level_start:level_end,n),surface_visits(n),local_key(4,n))
    visits=0;surface_visits=0;ncolumn=0;compile_level=0
    do d=1,size(grid)
       do p=3,grid(d)%patch%length
          call apply_onescale_to_patch(register_visit,grid(d),p-1,0,0,1)
       end do
    end do
    do l=level_end,level_start,-1
       compile_level=l
       do d=1,size(grid)
          do j=1,grid(d)%lev(l)%length
             call apply_onescale_to_patch(register_visit,grid(d),grid(d)%lev(l)%elts(j),0,0,1)
          end do
       end do
    end do
    allocate(sc(n_process),sd(n_process),rc(n_process),rd(n_process),counts(n_process))
    sc=0
    do c=1,ncolumn
       if(surface_visits(c)<1.or.sum(visits(:,c))<1) call die('inconsistent pressure visit closure')
       call inverse_operand_key(local_domain(c),local_node(c)-1,local_key(:,c))
       b=local_key(1,c)
       if(b==0) then
          ! An alias has a physical storage identity, not a canonical node.
          ! Assign this ordered column to the first catalogue block of its
          ! parent Domain. Never collapse aliases by coordinates.
          do b=1,size(block_catalog)
             if(block_catalog(b)%root_domain==glo_id(rank+1,local_domain(c))) exit
          end do
          if(b>size(block_catalog)) call die('alias has no parent block')
          local_key(:,c)=[b,0,0,0]
       end if
       r=block_catalog(b)%owner+1
       sc(r)=sc(r)+1
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr);call check(ierr,'column counts')
    call offsets(sc,sd);call offsets(rc,rd)
    h=5+level_end-level_start+1
    allocate(send_key(h,ncolumn),request(h,sum(rc)),order(ncolumn),packed_local(ncolumn))
    counts=0
    do c=1,ncolumn
       r=block_catalog(local_key(1,c))%owner+1;pos=sd(r)+counts(r)+1;counts(r)=counts(r)+1
       send_key(:,pos)=[local_key(:,c),surface_visits(c),visits(:,c)]
       order(c)=pos;packed_local(pos)=c
    end do
    call MPI_Alltoallv(send_key,h*sc,h*sd,MPI_INTEGER,request,h*rc,h*rd,MPI_INTEGER,comm,ierr)
    call check(ierr,'column address/visit plan')
    allocate(block_head(size(block_catalog)),next_column(sum(rc)))
    block_head=0;next_column=0
    do pos=size(request,2),1,-1
       if(request(2,pos)==0) cycle
       b=request(1,pos)
       if(block_catalog(b)%owner/=rank) call die('column delivered to wrong block owner')
       next_column(pos)=block_head(b);block_head(b)=pos
    end do
    allocate(alias_sc(n_process),alias_sd(n_process),alias_rc(n_process),alias_rd(n_process),alias_index(sum(rc)))
    alias_sc=0;alias_rc=0;alias_index=0
    do r=1,n_process
       do pos=sd(r)+1,sd(r)+sc(r)
          if(send_key(2,pos)==0) alias_sc(r)=alias_sc(r)+1
       end do
       do pos=rd(r)+1,rd(r)+rc(r)
          if(request(2,pos)==0) alias_rc(r)=alias_rc(r)+1
       end do
    end do
    call offsets(alias_sc,alias_sd);call offsets(alias_rc,alias_rd)
    a=0
    do pos=1,size(request,2)
       if(request(2,pos)/=0) cycle
       a=a+1;alias_index(pos)=a
    end do
    allocate(phase_sc(n_process),phase_sd(n_process),phase_rc(n_process),phase_rd(n_process), &
         first(ncolumn),cursor(ncolumn),remote_first(sum(rc)),lower(sum(rc)),geopot(sum(rc)),surface(ncolumn))
    generation=inverse_operand_generation()
    write(*,'(a,i0,a,4(1x,i0))') 'Block columns plan rank=',rank,' geometry/owner/aliases/domains=', &
         ncolumn,sum(rc),sum(alias_rc),size(grid)
  end subroutine

  subroutine start_block_columns
    integer :: c,p,d,id,k,a,ierr
    real(dp), allocatable :: send_alias(:,:,:),geo_send(:),geo_recv(:)
    type(Column_Context) :: context
    if(.not.armed) return
    call require_tendency_operands
    if(.not.compressible) call die('incompressible integration is not implemented')
    if(inverse_operand_generation()<0) call die('missing provisional inverse operands')
    if(generation/=inverse_operand_generation()) call compile_columns
    allocate(send_alias(4,zlevels,sum(alias_sc)),alias_operand(4,zlevels,sum(alias_rc)))
    allocate(geo_send(ncolumn),geo_recv(sum(rc)))
    a=0
    do p=1,ncolumn
       c=packed_local(p);d=local_domain(c);id=local_node(c)
       ! Surface geopotential is a case-specific reference input. It is not
       ! inferred from coordinates or silently replaced by zero topography.
       geo_send(p)=surf_geopot(d,id)
       if(local_key(2,c)/=0) cycle
       a=a+1
       do k=1,zlevels
          send_alias(:,k,a)=[inverse_scalar_operand(d,id-1,S_MASS,k), &
               inverse_scalar_operand(d,id-1,S_TEMP,k),sol_mean(S_MASS,k)%data(d)%elts(id), &
               sol_mean(S_TEMP,k)%data(d)%elts(id)]
       end do
    end do
    call MPI_Alltoallv(send_alias,4*zlevels*alias_sc,4*zlevels*alias_sd,MPI_DOUBLE_PRECISION, &
         alias_operand,4*zlevels*alias_rc,4*zlevels*alias_rd,MPI_DOUBLE_PRECISION,comm,ierr)
    call check(ierr,'inverse alias operands')
    call MPI_Alltoallv(geo_send,sc,sd,MPI_DOUBLE_PRECISION,geo_recv,rc,rd,MPI_DOUBLE_PRECISION,comm,ierr)
    call check(ierr,'surface reference operands')
    geopot=geo_recv
    allocate(filled(sum(rc)));filled=.false.
    context=Column_Context(0,0)
    call apply_local_block_field_consumer(compute_columns,context)
    call compute_aliases(context)
    if(.not.all(filled)) call die('missing final-owner column operands')
    deallocate(filled)
    call MPI_Alltoallv(lower,rc,rd,MPI_DOUBLE_PRECISION,surface,sc,sd,MPI_DOUBLE_PRECISION,comm,ierr)
    call check(ierr,'surface results')
    columns_active=.true.;armed=.false.;initial_pending=.false.;expected_layer=1;expected_level=level_end
    phase_layer=0;phase_level=0
  end subroutine

  subroutine compute_columns(catalog,block,context)
    integer, intent(in) :: catalog
    type(Block_Data), intent(in) :: block
    class(*), intent(inout) :: context
    integer :: p,node,n,k,mi,ti,sm,st,patch
    real(dp) :: operands(4,zlevels)
    select type(context)
    type is(Column_Context)
       n=size(block%node)
       if(block%scalar_mult/=1.or.block%field_level>1.or. &
            block%field_level+block%n_field_level-1<zlevels) &
            call die('unsupported resident scalar layout')
       sm=S_MASS-block%scalar_variable;st=S_TEMP-block%scalar_variable
       if(min(sm,st)<0.or.max(sm,st)>=block%n_scalar_variable) call die('missing resident variable')
       p=block_head(catalog)
       do while(p>0)
          if(request(2,p)/=STORE_PATCH) call die('unsupported compact column key')
          patch=request(3,p)
          if(patch<1.or.patch>size(block%patch)) call die('invalid column patch')
          node=block%patch(patch)%elts_start+request(4,p)+1
          if(node<1.or.node>n) call die('invalid compact column node')
          if(context%layer>0) then
             k=context%layer
             mi=(sm*block%n_field_level+k-block%field_level)*n+node;ti=(st*block%n_field_level+k-block%field_level)*n+node
             operands(:,k)=[block%scalar(mi),block%scalar(ti),block%scalar_mean(mi),block%scalar_mean(ti)]
          else
             do k=1,zlevels
                mi=(sm*block%n_field_level+k-block%field_level)*n+node;ti=(st*block%n_field_level+k-block%field_level)*n+node
                operands(:,k)=[block%scalar(mi),block%scalar(ti),block%scalar_mean(mi),block%scalar_mean(ti)]
             end do
          end if
          call evaluate_column(p,operands,context)
          p=next_column(p)
       end do
    class default
       call die('column callback context')
    end select
  end subroutine

  subroutine compute_aliases(context)
    type(Column_Context), intent(in) :: context
    integer :: p,a
    do p=1,size(request,2)
       a=alias_index(p)
       if(a==0) cycle
       call evaluate_column(p,alias_operand(:,:,a),context)
    end do
  end subroutine

  subroutine evaluate_column(p,operands,context)
    integer, intent(in) :: p
    real(dp), intent(in) :: operands(4,zlevels)
    type(Column_Context), intent(in) :: context
    integer :: k,j,pos,n
    real(dp) :: pressure,exner,geopot_lower
    if(filled(p)) call die('duplicate column producer')
    filled(p)=.true.
    k=context%layer
    if(k==0) then
       if(.not.all(ieee_is_finite(operands))) call die('nonfinite surface operand')
       lower(p)=0.0_dp
       do k=1,zlevels
          do j=1,request(5,p)
             lower(p)=lower(p)+(operands(3,k)+operands(1,k))
          end do
       end do
       lower(p)=grav_accel*lower(p)+p_top
    else
       if(.not.all(ieee_is_finite(operands(:,k)))) call die('nonfinite pressure operand')
       n=request(6+context%level-level_start,p)
       do j=1,n
          call pressure_layer(operands(1,k),operands(3,k),operands(2,k),operands(4,k), &
               grav_accel,c_p,p_0,kappa,lower(p),geopot(p),pressure,exner,geopot_lower)
          pos=remote_first(p)+j-1
          send_result(:,pos)=[lower(p),geopot(p),pressure,exner,geopot_lower]
          executed=executed+1_int64
       end do
    end if
  end subroutine

  subroutine column_consumer_address(d,id,owner_rank,address)
    integer, intent(in) :: d,id
    integer, intent(out) :: owner_rank,address
    integer :: c,r,b
    if(.not.columns_active) call die('column producer is not active')
    c=lookup(d)%column(id+1)
    if(c==0) then
       ! The regular traversal also emits unused Bernoulli boundary writes.
       ! Their PV/vorticity still execute on a parent block owner. Every
       ! Bernoulli consumer is guarded separately; no pressure fallback.
       do b=1,size(block_catalog)
          if(block_catalog(b)%root_domain==glo_id(rank+1,d)) exit
       end do
       if(b>size(block_catalog)) call die('boundary primitive has no parent block')
       owner_rank=block_catalog(b)%owner;address=0;return
    end if
    r=block_catalog(local_key(1,c))%owner+1
    owner_rank=r-1;address=order(c)-sd(r)
  end subroutine

  function owner_column_result(source_rank,address,visit,k,l) result(values)
    integer, intent(in) :: source_rank,address,visit,k,l
    real(dp) :: values(5)
    integer :: c,pos
    if(.not.columns_active.or.k/=phase_layer.or.l/=phase_level) call die('primitive pressure epoch')
    if(source_rank<0.or.source_rank>=n_process) call die('primitive pressure source')
    if(address<1.or.address>rc(source_rank+1)) call die('primitive pressure address')
    c=rd(source_rank+1)+address
    if(visit<1.or.visit>request(6+l-level_start,c)) call die('primitive pressure visit')
    pos=remote_first(c)+visit-1
    values=send_result(:,pos)
  end function

  subroutine column_surface(d,values)
    integer, intent(in) :: d
    real(dp), intent(inout) :: values(:)
    integer :: id,c,p
    if(.not.columns_active) return
    do id=1,size(lookup(d)%column)
       c=lookup(d)%column(id);if(c==0) cycle
       p=order(c)
       if(columns_oracle) call exact([values(id)],[surface(p)],'surface pressure')
       values(id)=surface(p)
    end do
  end subroutine

  subroutine finish_phase
    integer :: c
    if(phase_layer==0) return
    do c=1,ncolumn
       if(cursor(c)/=visits(phase_level,c)) call die('unconsumed pressure result')
    end do
  end subroutine

  subroutine begin_column_phase(k,l)
    integer, intent(in) :: k,l
    integer :: r,p,c,pos,ierr
    type(Column_Context) :: context
    if(.not.columns_active) return
    if(generation/=inverse_operand_generation()) call die('column topology changed during transaction')
    call finish_phase
    if(k/=expected_layer.or.l/=expected_level) call die('pressure producer phase order')
    expected_level=l-1
    if(l==level_start) then
       expected_layer=k+1;expected_level=level_end
    end if
    phase_layer=k;phase_level=l
    phase_sc=0;phase_rc=0
    do r=1,n_process
       do p=sd(r)+1,sd(r)+sc(r)
          c=packed_local(p);phase_sc(r)=phase_sc(r)+visits(l,c)
       end do
       do p=rd(r)+1,rd(r)+rc(r)
          phase_rc(r)=phase_rc(r)+request(6+l-level_start,p)
       end do
    end do
    call offsets(phase_sc,phase_sd);call offsets(phase_rc,phase_rd)
    if(allocated(send_result)) deallocate(send_result,result)
    allocate(send_result(5,sum(phase_rc)),result(5,sum(phase_sc)))
    pos=1
    do p=1,size(request,2)
       remote_first(p)=pos;pos=pos+request(6+l-level_start,p)
    end do
    pos=1
    do p=1,ncolumn
       c=packed_local(p);first(c)=pos;pos=pos+visits(l,c)
    end do
    cursor=0
    allocate(filled(sum(rc)));filled=.false.
    context=Column_Context(k,l)
    call apply_local_block_field_consumer(compute_columns,context)
    call compute_aliases(context)
    if(.not.all(filled)) call die('missing pressure producer')
    deallocate(filled)
    call MPI_Alltoallv(send_result,5*phase_rc,5*phase_rd,MPI_DOUBLE_PRECISION, &
         result,5*phase_sc,5*phase_sd,MPI_DOUBLE_PRECISION,comm,ierr)
    call check(ierr,'pressure result delivery')
  end subroutine

  subroutine exact(actual,expected,where)
    real(dp), intent(in) :: actual(:),expected(:)
    character(*), intent(in) :: where
    if(size(actual)/=size(expected)) call die('oracle extent')
    if(any(transfer(actual,[0_int64],size(actual))/=transfer(expected,[0_int64],size(expected)))) then
       write(*,*) 'actual:',actual,' expected:',expected
       call die('independent numerical mismatch: '//where)
    end if
    checked=checked+1_int64
  end subroutine

  subroutine column_result(d,id,k,pl,g,p,e,gl)
    integer, intent(in) :: d,id,k
    real(dp), intent(inout) :: pl,g,p,e,gl
    integer :: c,pos
    if(.not.columns_active) call die('pressure consumer outside transaction')
    if(k/=phase_layer) call die('pressure consumer layer')
    c=lookup(d)%column(id)
    if(c==0) call die('pressure consumer is missing from plan')
    cursor(c)=cursor(c)+1
    if(cursor(c)>visits(phase_level,c)) call die('extra pressure consumer visit')
    pos=first(c)+cursor(c)-1
    if(columns_oracle) call exact([pl,g,p,e,gl],result(:,pos),'pressure layer')
    pl=result(1,pos);g=result(2,pos);p=result(3,pos);e=result(4,pos);gl=result(5,pos)
  end subroutine

  subroutine end_block_columns
    if(.not.columns_active) return
    call finish_phase
    if(expected_layer/=zlevels+1.or.expected_level/=level_end) call die('incomplete column transaction')
    deallocate(alias_operand,send_result,result)
    columns_active=.false.;phase_layer=0
    write(*,'(a,i0,a,2(1x,i0))') 'Block columns complete rank=',rank,' executed/checked=',executed,checked
  end subroutine
end module
