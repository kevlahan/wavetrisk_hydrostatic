module parallel_block_rhs_retention_mod
  ! The authoritative connected RHS is retained directly in the existing
  ! resident RK accumulator. Only topology identities cross the geometry cut.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : zlevels,scalars,NORTH,EAST,NORTHEAST,p_top,grav_accel,c_p,p_0,kappa
  use arch_mod, only : rank,glo_id
  use domain_mod, only : grid
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data,STORE_PATCH,STORE_BDRY,STORE_GHOST,apply_local_block_field_consumer, &
       apply_local_block_tendency_kernel,retain_local_block_tendency_in_accumulator, &
       n_local_blocks,local_block_catalog,local_block_patch_count,catalog_local_block
  use parallel_block_rhs_placement_mod, only : compile_rhs_owner_map
  use parallel_block_rhs_distribute_mod, only : deliver_rhs_metadata
  use parallel_block_rhs_distributed_mod, only : Rhs_Distributed_Executor
  implicit none
  private
  public :: Rhs_Retention,prepare_rhs_retention,retain_rhs_tendency
  type :: Rhs_Rounding_Column
     logical,allocatable :: ready(:)
     real(dp),allocatable :: exner(:,:),theta(:,:)
  end type
  type :: Block_Rhs_Index
     integer :: catalog=0
     integer,allocatable :: node(:)
     real(dp),allocatable :: length(:,:)
     type(Rhs_Rounding_Column) :: column(3)
  end type
  type :: Rhs_Retention
     type(Block_Rhs_Index),allocatable :: block(:)
     type(Rhs_Distributed_Executor),pointer :: executor=>null()
     integer(int64) :: retained_values=0,retained_stages=0
  end type
contains
  subroutine prepare_rhs_retention(executor,base,topology,retention)
    type(Rhs_Distributed_Executor),intent(in) :: executor
    integer,intent(in) :: base(:)
    integer(int64),intent(in) :: topology
    type(Rhs_Retention),intent(out) :: retention
    integer,allocatable :: physical(:,:),state(:,:),owner(:),key(:,:)
    real(dp),allocatable :: resolved(:,:),lengths(:,:),physical_lengths(:,:)
    integer :: d,b,n,i,c,p,j,local
    n=base(size(base));allocate(physical(4,n),state(4,n))
    do d=1,size(grid)
       b=base(d)
       call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),owner,key)
       do i=1,grid(d)%node%length
          physical(:,b+i)=[glo_id(rank+1,d)+1,i,1,1]
       end do
       state(:,b+1:b+size(key,2))=key
    end do
    call deliver_rhs_metadata(physical,real(transpose(state),dp),executor%placement%key,resolved,topology,8401)
    if(any(abs(real(nint(resolved),dp)-resolved)>0.0_dp)) error stop 'RHS retention identity'
    state=transpose(nint(resolved))
    allocate(physical_lengths(n,3))
    do d=1,size(grid)
       do i=1,grid(d)%node%length
          physical_lengths(base(d)+i,:)=grid(d)%len%elts(3*i-2:3*i)
       end do
    end do
    call deliver_rhs_metadata(physical,physical_lengths,executor%placement%key,resolved,topology,8402)
    allocate(lengths(size(resolved,2),size(resolved,1)))
    lengths(:,:)=transpose(resolved)
    allocate(retention%block(n_local_blocks()))
    do local=1,size(retention%block)
       c=local_block_catalog(local);retention%block(local)%catalog=c
       allocate(retention%block(local)%node(PATCH_SIZE**2*local_block_patch_count(c)))
       retention%block(local)%node=0
    end do
    do i=1,size(state,2)
       if(.not.executor%placement%owned(i)) cycle
       c=state(1,i)
       if(c==0) cycle
       if(state(2,i)/=STORE_PATCH) error stop 'RHS retention noninterior key'
       local=catalog_local_block(c)
       if(local<1.or.local>size(retention%block)) error stop 'RHS retention wrong owner'
       p=state(3,i);j=PATCH_SIZE**2*(p-1)+state(4,i)+1
       if(j<1.or.j>size(retention%block(local)%node)) error stop 'RHS retention node extent'
       if(retention%block(local)%node(j)/=0) error stop 'RHS retention duplicate interior'
       retention%block(local)%node(j)=i
    end do
    do local=1,size(retention%block)
       if(any(retention%block(local)%node==0)) error stop 'RHS retention incomplete interior'
       retention%block(local)%length=lengths(:,retention%block(local)%node)
       if(any(retention%block(local)%length<=0.0_dp)) error stop 'RHS retention missing edge metric'
    end do
    call apply_local_block_field_consumer(allocate_columns,retention)
  end subroutine

  subroutine retain_rhs_tendency(executor,retention)
    type(Rhs_Distributed_Executor),intent(in),target :: executor
    type(Rhs_Retention),intent(inout) :: retention
    retention%executor=>executor
    call apply_local_block_tendency_kernel(install,retention)
    call retain_local_block_tendency_in_accumulator
    nullify(retention%executor)
    retention%retained_stages=retention%retained_stages+1
  end subroutine

  subroutine install(catalog,block,scalar_tendency,vector_tendency,context)
    integer,intent(in) :: catalog
    type(Block_Data),intent(in) :: block
    real(dp),intent(inout) :: scalar_tendency(:),vector_tendency(:)
    class(*),intent(inout) :: context
    integer :: local,n,p,j,q,k,v,lev,a,b
    select type(c=>context)
    type is(Rhs_Retention)
       if(.not.associated(c%executor)) error stop 'RHS retention missing completed executor'
       local=catalog_local_block(catalog);n=size(block%node)
       if(c%block(local)%catalog/=catalog.or.size(c%block(local)%node)/=PATCH_SIZE**2*size(block%patch)) &
            error stop 'RHS retention stale block layout'
       if(block%scalar_variable/=scalars(1).or.block%n_scalar_variable/=2.or.block%scalar_mult/=1.or. &
            block%vector_mult/=3) error stop 'RHS retention unsupported family'
       do p=1,size(block%patch)
          do j=1,PATCH_SIZE**2
             q=block%patch(p)%elts_start+j
             a=c%block(local)%node(PATCH_SIZE**2*(p-1)+j)
             do k=1,zlevels
                lev=k-block%field_level
                if(lev<0.or.lev>=block%n_field_level) error stop 'RHS retention layer extent'
                do v=1,2
                   b=((v-1)*block%n_field_level+lev)*n+q
                   scalar_tendency(b)=c%executor%scalar_rhs(a,v,k)
                end do
                b=3*(lev*n+q-1)
                vector_tendency(b+1:b+3)=c%executor%velocity_rhs(3*a-2:3*a,k)
             end do
          end do
       end do
       call accepted_velocity_rounding(block,c%block(local),vector_tendency)
       c%retained_values=c%retained_values+int(5*PATCH_SIZE**2*size(block%patch),int64)*zlevels
    class default
       error stop 'RHS retention context'
    end select
  end subroutine
  subroutine allocate_columns(catalog,block,context)
    integer,intent(in) :: catalog
    type(Block_Data),intent(in) :: block
    class(*),intent(inout) :: context
    integer :: local,s,n
    select type(c=>context)
    type is(Rhs_Retention)
       local=catalog_local_block(catalog)
       do s=1,3
          select case(s)
          case(STORE_PATCH)
             n=size(block%node)
          case(STORE_BDRY)
             n=size(block%bdry_node)
          case(STORE_GHOST)
             n=size(block%ghost_node)
          end select
          allocate(c%block(local)%column(s)%ready(n),c%block(local)%column(s)%exner(zlevels,n), &
               c%block(local)%column(s)%theta(zlevels,n))
          c%block(local)%column(s)%ready=.false.
       end do
    class default
       error stop 'RHS rounding column context'
    end select
  end subroutine

  subroutine accepted_velocity_rounding(block,index,tendency)
    ! Accepted266 retains (raw_rhs + reconstructed_Exner_term) - term.
    ! The cancellation changes rounding. Preserve it on resident owners until
    ! a separately approved numerical baseline permits its removal.
    type(Block_Data),intent(in) :: block
    type(Block_Rhs_Index),intent(inout) :: index
    real(dp),intent(inout) :: tendency(:)
    integer :: p,i,j,q,k,e,node,neighbor,storage,s,b,link,record,dims(2),start
    real(dp) :: difference,remainder,lower,pressure,rho(zlevels),temperature(zlevels)
    do s=1,3
       index%column(s)%ready=.false.
    end do
    do p=1,size(block%patch)
       do q=0,PATCH_SIZE**2-1
          i=mod(q,PATCH_SIZE);j=q/PATCH_SIZE;node=block%patch(p)%elts_start+q
          call column(STORE_PATCH,node)
          do e=1,3
             storage=STORE_PATCH;neighbor=node
             select case(e)
             case(1)
                if(i<PATCH_SIZE-1) then
                   neighbor=node+1
                else
                   call adjacent(EAST,0,j)
                end if
             case(2)
                if(i<PATCH_SIZE-1.and.j<PATCH_SIZE-1) then
                   neighbor=node+PATCH_SIZE+1
                else if(i==PATCH_SIZE-1.and.j<PATCH_SIZE-1) then
                   call adjacent(EAST,0,j+1)
                else if(i<PATCH_SIZE-1.and.j==PATCH_SIZE-1) then
                   call adjacent(NORTH,i+1,0)
                else
                   call adjacent(NORTHEAST,0,0)
                end if
             case(3)
                if(j<PATCH_SIZE-1) then
                   neighbor=node+PATCH_SIZE
                else
                   call adjacent(NORTH,i,0)
                end if
             end select
             call column(storage,neighbor)
             do k=1,zlevels
                difference=index%column(storage)%exner(k,neighbor+1)-index%column(STORE_PATCH)%exner(k,node+1)
                if(e==2) difference=index%column(STORE_PATCH)%exner(k,node+1)-index%column(storage)%exner(k,neighbor+1)
                difference=difference*(0.5_dp*(index%column(STORE_PATCH)%theta(k,node+1)+ &
                     index%column(storage)%theta(k,neighbor+1)))
                difference=difference/index%length(e,PATCH_SIZE**2*(p-1)+q+1)
                b=3*((k-block%field_level)*size(block%node)+node)+e
                remainder=tendency(b)+difference
                tendency(b)=remainder-difference
             end do
          end do
       end do
    end do
  contains
    subroutine adjacent(requested_side,requested_i,requested_j)
      integer,intent(in) :: requested_side,requested_i,requested_j
      integer :: n
      n=block%patch(p)%neigh(requested_side)
      if(n>0) then
         storage=STORE_PATCH;start=block%patch(n+1)%elts_start;dims=PATCH_SIZE
      else if(n<0) then
         link=-n
         if(block%block_bdry(link)%patch/=p-1.or.block%block_bdry(link)%side/=requested_side) &
              error stop 'RHS rounding neighbor link'
         record=block%block_bdry(link)%ghost_id
         if(record>0) then
            storage=STORE_GHOST;start=block%ghost_storage(record)%local_start;dims=PATCH_SIZE
         else
            record=block%block_bdry(link)%storage_id
            storage=STORE_BDRY;start=block%bdry_storage(record)%local_start;dims=block%block_bdry(link)%dims
         end if
      else
         error stop 'RHS rounding absent neighbor'
      end if
      if(requested_i<0.or.requested_i>=dims(1).or.requested_j<0.or.requested_j>=dims(2)) &
           error stop 'RHS rounding neighbor coordinates'
      neighbor=start+requested_j*dims(1)+requested_i
    end subroutine
    subroutine column(storage,node)
      integer,intent(in) :: storage,node
      integer :: k
      if(index%column(storage)%ready(node+1)) return
      do k=1,zlevels
         rho(k)=total(storage,node,k,0);temperature(k)=total(storage,node,k,1)
      end do
      lower=p_top
      do k=1,zlevels
         if(rho(k)<=0.0_dp) error stop 'RHS rounding column mass'
         lower=lower+grav_accel*rho(k)
      end do
      do k=1,zlevels
         pressure=lower-0.5_dp*grav_accel*rho(k)
         if(pressure<=0.0_dp) error stop 'RHS rounding column pressure'
         index%column(storage)%exner(k,node+1)=c_p*(pressure/p_0)**kappa
         index%column(storage)%theta(k,node+1)=temperature(k)/rho(k)
         if(k<zlevels) lower=lower-grav_accel*rho(k)
      end do
      index%column(storage)%ready(node+1)=.true.
    end subroutine
    real(dp) function total(storage,node,k,v) result(value)
      integer,intent(in) :: storage,node,k,v
      integer :: b,n
      select case(storage)
      case(STORE_PATCH)
         n=size(block%node);b=(v*block%n_field_level+k-block%field_level)*n+node+1
         value=block%scalar(b)+block%scalar_mean(b)
      case(STORE_BDRY)
         n=size(block%bdry_node);b=(v*block%n_field_level+k-block%field_level)*n+node+1
         value=block%bdry_scalar(b)+block%bdry_scalar_mean(b)
      case(STORE_GHOST)
         n=size(block%ghost_node);b=(v*block%n_field_level+k-block%field_level)*n+node+1
         value=block%ghost_scalar(b)+block%ghost_scalar_mean(b)
      case default
         error stop 'RHS rounding column storage'
      end select
    end function
  end subroutine
end module
