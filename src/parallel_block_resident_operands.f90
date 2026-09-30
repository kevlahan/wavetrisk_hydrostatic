module parallel_block_resident_operands_mod
  ! Explicit stencil inputs: interiors come from final-owner Block_Data;
  ! physical aliases come from the named completed inverse epoch. No q argument.
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use ieee_arithmetic, only: ieee_is_finite
  use arch_mod, only: rank,n_process,comm,block_catalog,abort_run
  use shared_mod, only: S_MASS,S_TEMP,S_VELO,zlevels,EDGE
  use domain_mod, only: sol_mean
  use parallel_block_mod, only: Block_Data,STORE_PATCH,apply_local_block_field_consumer
  use parallel_block_inverse_mod, only: inverse_operand_key,inverse_scalar_operand,inverse_vector_operand, &
       require_tendency_operands
  implicit none
  private
  public :: Resident_Operand_Plan,compile_resident_operands,fetch_resident_operands,resident_operand_key
  type :: Resident_Operand_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),source(:,:),address(:),head(:),next(:)
     real(dp), allocatable :: send(:,:),recv(:,:)
     logical, allocatable :: ready(:)
     integer :: layer=0
  end type
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': resident operands: ',message
    call abort_run
    error stop 'resident operands abort returned'
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
  function resident_operand_key(d,id) result(key)
    integer, intent(in) :: d,id
    integer :: key(5),compact(4)
    call inverse_operand_key(d,id,compact)
    if(compact(1)>0) then
       key=[block_catalog(compact(1))%owner,compact]
    else
       key=[rank,0,d,id,0]
    end if
  end function
  subroutine compile_resident_operands(plan,keys)
    type(Resident_Operand_Plan), intent(out) :: plan
    integer, intent(in) :: keys(:,:)
    integer :: n,p,r,pos,ierr,b,nunique,t,h,slot,u
    integer(int64) :: hash
    integer, allocatable :: packed(:,:),count(:),unique(:,:),table(:),mapping(:),packed_unique(:)
    if(size(keys,1)/=5) call die('key extent')
    n=size(keys,2)
    ! One transfer per distinct physical operand, even when many operations
    ! share the stencil. Aliases keep their storage identity in all five keys.
    allocate(unique(5,n),mapping(n),table(max(1,2*n+1)),packed_unique(n))
    table=0;nunique=0
    do p=1,n
       hash=0_int64
       do t=1,5
          hash=modulo(65599_int64*hash+int(keys(t,p),int64),2147483647_int64)
       end do
       h=int(modulo(hash,int(size(table),int64)))+1
       do
          slot=table(h)
          if(slot==0) then
             nunique=nunique+1;unique(:,nunique)=keys(:,p);table(h)=nunique;slot=nunique;exit
          end if
          if(all(unique(:,slot)==keys(:,p))) exit
          h=modulo(h,size(table))+1
       end do
       mapping(p)=slot
    end do
    allocate(plan%sc(n_process),plan%sd(n_process),plan%rc(n_process),plan%rd(n_process),count(n_process))
    plan%sc=0
    do p=1,nunique
       r=unique(1,p)+1
       if(r<1.or.r>n_process) call die('source rank')
       plan%sc(r)=plan%sc(r)+1
    end do
    call MPI_Alltoall(plan%sc,1,MPI_INTEGER,plan%rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('request counts')
    call offsets(plan%sc,plan%sd);call offsets(plan%rc,plan%rd)
    allocate(packed(4,nunique),plan%source(4,sum(plan%rc)),plan%address(n))
    count=0
    do p=1,nunique
       r=unique(1,p)+1;pos=plan%sd(r)+count(r)+1;count(r)=count(r)+1
       packed(:,pos)=unique(2:5,p);packed_unique(p)=pos
    end do
    do p=1,n
       u=mapping(p);plan%address(p)=packed_unique(u)
    end do
    call MPI_Alltoallv(packed,4*plan%sc,4*plan%sd,MPI_INTEGER,plan%source,4*plan%rc,4*plan%rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('operand requests')
    allocate(plan%head(size(block_catalog)),plan%next(sum(plan%rc)))
    plan%head=0;plan%next=0
    do p=size(plan%source,2),1,-1
       b=plan%source(1,p)
       if(b==0) cycle
       if(b<1.or.b>size(block_catalog)) call die('catalogue index')
       if(block_catalog(b)%owner/=rank) call die('operand request sent to wrong owner')
       if(plan%source(2,p)/=STORE_PATCH) call die('unsupported resident storage key')
       plan%next(p)=plan%head(b);plan%head(b)=p
    end do
    allocate(plan%send(7,sum(plan%rc)),plan%recv(7,nunique),plan%ready(sum(plan%rc)))
  end subroutine
  subroutine fetch_resident_operands(plan,k)
    type(Resident_Operand_Plan), intent(inout) :: plan
    integer, intent(in) :: k
    integer :: p,d,id,ierr
    call require_tendency_operands
    if(k<1.or.k>zlevels) call die('layer')
    plan%layer=k;plan%ready=.false.
    call apply_local_block_field_consumer(read_block,plan)
    do p=1,size(plan%source,2)
       if(plan%source(1,p)/=0) cycle
       d=plan%source(2,p);id=plan%source(3,p)
       plan%send(1,p)=inverse_scalar_operand(d,id,S_MASS,k)+sol_mean(S_MASS,k)%data(d)%elts(id+1)
       plan%send(2,p)=inverse_scalar_operand(d,id,S_TEMP,k)+sol_mean(S_TEMP,k)%data(d)%elts(id+1)
       plan%send(3:5,p)=inverse_vector_operand(d,id,k)
       plan%send(6,p)=inverse_scalar_operand(d,id,S_MASS,k)
       plan%send(7,p)=inverse_scalar_operand(d,id,S_TEMP,k)
       plan%ready(p)=.true.
    end do
    if(.not.all(plan%ready)) call die('missing resident operand')
    if(.not.all(ieee_is_finite(plan%send))) call die('nonfinite resident operand')
    call MPI_Alltoallv(plan%send,7*plan%rc,7*plan%rd,MPI_DOUBLE_PRECISION, &
         plan%recv,7*plan%sc,7*plan%sd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('operand delivery')
  end subroutine
  subroutine read_block(catalog,block,context)
    integer, intent(in) :: catalog
    type(Block_Data), intent(in) :: block
    class(*), intent(inout) :: context
    integer :: p,patch,node,n,k,mi,ti,vi,sm,st
    select type(context)
    type is(Resident_Operand_Plan)
       p=context%head(catalog)
       if(p==0) return
       k=context%layer;n=size(block%node)
       if(block%scalar_mult/=1.or.block%vector_mult/=EDGE.or.block%vector_variable/=S_VELO) &
            call die('resident field layout')
       if(k<block%field_level.or.k>=block%field_level+block%n_field_level) call die('resident layer coverage')
       sm=S_MASS-block%scalar_variable;st=S_TEMP-block%scalar_variable
       if(min(sm,st)<0.or.max(sm,st)>=block%n_scalar_variable) call die('resident scalar coverage')
       do while(p>0)
          patch=context%source(3,p)
          if(patch<1.or.patch>size(block%patch)) call die('resident patch')
          node=block%patch(patch)%elts_start+context%source(4,p)+1
          if(node<1.or.node>n) call die('resident node')
          mi=(sm*block%n_field_level+k-block%field_level)*n+node
          ti=(st*block%n_field_level+k-block%field_level)*n+node
          vi=EDGE*((k-block%field_level)*n+node-1)
          context%send(1,p)=block%scalar(mi)+block%scalar_mean(mi)
          context%send(2,p)=block%scalar(ti)+block%scalar_mean(ti)
          context%send(3:5,p)=block%vector(vi+1:vi+EDGE)
          context%send(6,p)=block%scalar(mi);context%send(7,p)=block%scalar(ti)
          context%ready(p)=.true.;p=context%next(p)
       end do
    class default
       call die('resident operand callback')
    end select
  end subroutine
end module
