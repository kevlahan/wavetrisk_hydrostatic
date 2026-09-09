module parallel_block_inverse_mod
  ! Native inverse boundary dependency workspace. Topology is compiled once
  ! per generation; only required nodes are staged, never Float_Field copies.
  ! Regular scalar/inner-vector kernels remain on final-owner blocks. The
  ! sparse outer-edge producer retains the established geometric route owner
  ! and write order, including aliases which have no compact interior owner.
  use iso_fortran_env, only : int64, error_unit
  use ieee_arithmetic, only : ieee_is_finite
  use mpi_f08, only : MPI_Alltoall, MPI_Alltoallv, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUCCESS, MPI_Wtime, &
       MPI_Request, MPI_Irecv, MPI_Isend, MPI_Waitall, MPI_STATUSES_IGNORE
  use kind_mod, only : dp
  use shared_mod, only : EDGE, AT_NODE, AT_EDGE, N_VARIABLE, N_BDRY, N_CHDRN, &
       scalars, S_VELO, zlevels, n_domain, level_start, level_end, end_pt, opp_no, hex_sides, hex_s_offs, &
       RT, UP, NORTHEAST, NORTHWEST, IMINUSJPLUS, IPLUSJMINUS, IJMINUS
  use domain_mod, only : grid, sol, Float_Field, idx, ed_idx, get_offs_Domain, chd_offs
  use patch_mod, only : PATCH_SIZE, LAST
  use arch_mod, only : comm, rank, n_process, glo_id, loc_id, owner, block_catalog, abort_run
  use parallel_block_mod, only : Block_Data, STORE_PATCH, STORE_BDRY, BLOCK_PAYLOAD_SOL, &
       BLOCK_PAYLOAD_WAV_COEFF, apply_local_block_field_consumer, transfer_local_block_inverse_node
  implicit none
  private
  public :: prepare_native_inverse, native_inverse_gather, native_inverse_boundary, &
       native_inverse_scatter, native_inverse_outer, compare_native_inverse, publish_native_inverse
  real(dp), public :: native_inverse_seconds(3)=0.0_dp
  integer(int64), public :: native_inverse_calls(3)=0_int64, native_inverse_messages(3)=0_int64, &
       native_inverse_bytes(3)=0_int64
  logical :: timing=.false.

  type :: Node_Map
     integer, allocatable :: slot(:), key(:,:)
  end type
  type :: Transfer_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),key(:,:),slot(:)
     real(dp), allocatable :: send(:),recv(:)
  end type
  type :: Alias_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),source(:),dest(:)
     real(dp), allocatable :: send(:),recv(:)
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
  type(Alias_Plan), allocatable :: aliases(:,:)
  type(Local_Aliases) :: local_alias(2)
  type(Outer_Operation), allocatable :: operations(:)
  type(MPI_Request), allocatable :: requests(:)
  integer, allocatable :: node_domain(:),node_id(:),operation_first(:),operation_last(:)
  integer, allocatable :: scaffold_slot(:)
  real(dp), allocatable :: scaffold_value(:,:,:)
  integer(int64), allocatable :: coverage(:,:)
  real(dp), allocatable :: value(:,:,:,:)
  integer(int64) :: plan_generation=-1_int64
  integer :: nnode=0,noperation=0,nscalar=0
  real(dp), parameter :: base_weight(9) = &
       [16.0_dp,-1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,-1.0_dp,1.0_dp,1.0_dp]/16.0_dp

contains

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

  subroutine prepare_native_inverse(generation,scaling,wavelet,first_level,profile)
    integer(int64), intent(in) :: generation
    integer, intent(in) :: first_level
    logical, intent(in) :: profile
    type(Float_Field), intent(in) :: scaling(1:N_VARIABLE,1:zlevels),wavelet(1:N_VARIABLE,1:zlevels)
    integer :: d,id,q,k,v,total,b,p,r,pos,n,next_patch
    integer, allocatable :: sc(:),rc(:),sd(:),rd(:),sb(:,:),rb(:,:),request(:,:),slots(:)
    type(Boundary_Manifest) :: manifest

    timing=profile
    native_inverse_seconds=0.0_dp
    native_inverse_calls=0_int64
    native_inverse_messages=0_int64
    native_inverse_bytes=0_int64
    if (plan_generation /= generation) then
       if (allocated(nodes)) deallocate(nodes,node_domain,node_id,value,aliases,operations,coverage, &
            scaffold_slot,scaffold_value,operation_first,operation_last)
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
       allocate(aliases(2,level_start-1:level_end))
       call compile_aliases
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
       allocate(value(nscalar+EDGE,zlevels,nnode,2))
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(2,node_id(p)+1)==-1) n=n+1
       end do
       allocate(scaffold_slot(n),scaffold_value(nscalar+EDGE,zlevels,n))
       n=0
       do p=1,nnode
          if (nodes(node_domain(p))%key(2,node_id(p)+1)/=-1) cycle
          n=n+1
          scaffold_slot(n)=p
       end do
       plan_generation=generation
    end if

    ! Seed only compatibility/scaffold and pre-existing boundary state once
    ! per transform. Every catalogued interior dependency is replaced below
    ! from the authoritative final owner before any arithmetic or alias copy.
    do p=1,nnode
       d=node_domain(p)
       id=node_id(p)
       if (nodes(d)%key(1,id+1)>0) cycle
       do k=1,zlevels
          do v=1,nscalar
             value(v,k,p,1)=scaling(scalars(1)+v-1,k)%data(d)%elts(id+1)
             value(v,k,p,2)=wavelet(scalars(1)+v-1,k)%data(d)%elts(id+1)
          end do
          do v=1,EDGE
             value(nscalar+v,k,p,1)=scaling(S_VELO,k)%data(d)%elts(EDGE*id+v)
             value(nscalar+v,k,p,2)=wavelet(S_VELO,k)%data(d)%elts(EDGE*id+v)
          end do
       end do
    end do
    ! Legacy writeback preserves uncatalogued interior patches from global
    ! sol, NOT from the provisional scaling argument. Keep this fixed input
    ! separate from the pre-boundary seed and apply it at each gather point.
    do q=1,size(scaffold_slot)
       p=scaffold_slot(q)
       d=node_domain(p)
       id=node_id(p)
       do k=1,zlevels
          do v=1,nscalar
             scaffold_value(v,k,q)=sol(scalars(1)+v-1,k)%data(d)%elts(id+1)
          end do
          do v=1,EDGE
             scaffold_value(nscalar+v,k,q)=sol(S_VELO,k)%data(d)%elts(EDGE*id+v)
          end do
       end do
    end do
    do q=1,2
       do b=1,2
          call native_inverse_gather(q,b,.false.)
       end do
    end do
    do q=1,2
       call native_inverse_boundary(q,2,level_start,level_end)
       call native_inverse_boundary(q,1,first_level,level_end)
       call native_inverse_scatter(q,2,.false.)
       call native_inverse_scatter(q,1,.false.)
    end do

  contains
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
    integer, allocatable :: send_key(:,:)
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
    allocate(send_key(4,max(1,sum(plan%sc))),plan%key(4,max(1,sum(plan%rc))))
    allocate(plan%slot(sum(plan%sc)))
    plan%sc=0
    do p=1,size(slots)
       r=block_catalog(request(1,p))%owner+1
       pos=plan%sd(r)+plan%sc(r)+1
       send_key(:,pos)=request(:,p)
       plan%slot(pos)=slots(p)
       plan%sc(r)=plan%sc(r)+1
    end do
    call MPI_Alltoallv(send_key,4*plan%sc,4*plan%sd,MPI_INTEGER, &
         plan%key,4*plan%rc,4*plan%rd,MPI_INTEGER,comm,ierr)
    call check(ierr,'native dependency keys')
    allocate(plan%send(max(1,max(nscalar,EDGE)*zlevels*sum(plan%sc))))
    allocate(plan%recv(max(1,max(nscalar,EDGE)*zlevels*sum(plan%rc))))
  end subroutine

  subroutine native_inverse_gather(component,family,reset_scaffold)
    integer, intent(in) :: component,family
    logical, optional, intent(in) :: reset_scaffold
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
    call transfer_nodes(interior_plan,component,family,.false.)
  end subroutine

  subroutine native_inverse_scatter(component,family,interiors)
    integer, intent(in) :: component,family
    logical, intent(in) :: interiors
    if (interiors) call transfer_nodes(interior_plan,component,family,.true.)
    call transfer_nodes(boundary_plan,component,family,.true.)
  end subroutine

  subroutine transfer_nodes(plan,component,family,install)
    type(Transfer_Plan), intent(inout) :: plan
    integer, intent(in) :: component,family
    logical, intent(in) :: install
    integer :: p,q,k,v,nv,first,n,phase
    real(dp) :: started
    real(dp) :: sample(merge(nscalar,EDGE,component==1),zlevels)
    first=merge(1,nscalar+1,component==1)
    nv=size(sample,1)
    n=nv*zlevels
    started=0.0_dp
    if (timing) started=MPI_Wtime()
    if (install) then
       do p=1,sum(plan%sc)
          q=(p-1)*n
          do k=1,zlevels
             do v=1,nv
                q=q+1
                plan%send(q)=value(first+v-1,k,plan%slot(p),family)
             end do
          end do
       end do
       call exchange_values(plan%send,plan%sc,plan%sd,plan%recv,plan%rc,plan%rd,n,19071)
       do p=1,sum(plan%rc)
          sample=reshape(plan%recv((p-1)*n+1:p*n),shape(sample))
          call transfer_local_block_inverse_node(plan%key(:,p),family,component,.true.,sample)
       end do
    else
       do p=1,sum(plan%rc)
          call transfer_local_block_inverse_node(plan%key(:,p),family,component,.false.,sample)
          plan%recv((p-1)*n+1:p*n)=reshape(sample,[n])
       end do
       call exchange_values(plan%recv,plan%rc,plan%rd,plan%send,plan%sc,plan%sd,n,19072)
       do p=1,sum(plan%sc)
          q=(p-1)*n
          do k=1,zlevels
             do v=1,nv
                q=q+1
                value(first+v-1,k,plan%slot(p),family)=plan%send(q)
             end do
          end do
       end do
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

  subroutine exchange_values(send,sc,sd,recv,rc,rd,n,tag)
    ! Persistent sparse peer schedules, with no per-phase global collective.
    real(dp), intent(in), asynchronous :: send(:)
    real(dp), intent(out), asynchronous :: recv(:)
    integer, intent(in) :: sc(:),sd(:),rc(:),rd(:),n,tag
    integer :: r,nrequest,ierr
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
    recv(n*rd(r)+1:n*(rd(r)+rc(r)))=send(n*sd(r)+1:n*(sd(r)+sc(r)))
    if (nrequest>0) then
       call MPI_Waitall(nrequest,requests(1:nrequest),MPI_STATUSES_IGNORE,ierr)
       call check(ierr,'native dependency completion')
    end if
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
            allocate(a%send(max(1,nscalar*zlevels*sum(a%sc))),a%recv(max(1,nscalar*zlevels*sum(a%rc))))
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

  subroutine native_inverse_boundary(component,family,first,last)
    integer, intent(in) :: component,family,first,last
    integer :: l,p,k,v,nv,n,q
    real(dp) :: started
    started=0.0_dp
    if (timing) started=MPI_Wtime()
    nv=merge(nscalar,1,component==1)
    n=nv*zlevels
    ! Remote sends are snapshots BEFORE the ordered, all-level local copies.
    do l=first,last
       associate(a=>aliases(component,l))
         q=0
         do p=1,size(a%source)
            do k=1,zlevels
               do v=1,nv
                  q=q+1
                  a%send(q)=alias_value(a%source(p),component,v,k,family)
               end do
            end do
         end do
       end associate
    end do
    do p=1,size(local_alias(component)%source)
       do k=1,zlevels
          do v=1,nv
             call set_alias(local_alias(component)%dest(p),component,v,k,family, &
                  alias_value(local_alias(component)%source(p),component,v,k,family))
          end do
       end do
    end do
    do l=first,last
       associate(a=>aliases(component,l))
         call exchange_values(a%send,a%sc,a%sd,a%recv,a%rc,a%rd,n,19073)
         if (timing) call record_traffic(2,a%sc,n)
         q=0
         do p=1,size(a%dest)
            do k=1,zlevels
               do v=1,nv
                  q=q+1
                  call set_alias(a%dest(p),component,v,k,family,a%recv(q))
               end do
            end do
         end do
       end associate
    end do
    if (timing) then
       native_inverse_seconds(2)=native_inverse_seconds(2)+MPI_Wtime()-started
       native_inverse_calls(2)=native_inverse_calls(2)+1_int64
    end if
  end subroutine

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

  subroutine compare_native_inverse(scaling,component,payload_family)
    type(Float_Field), intent(in) :: scaling(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: component
    integer, optional, intent(in) :: payload_family
    integer :: p,d,id,k,v,first,last,var,address,family
    real(dp) :: reference,native,allowed
    first=merge(1,nscalar+1,component==1)
    last=merge(nscalar,nscalar+EDGE,component==1)
    family=BLOCK_PAYLOAD_SOL
    if (present(payload_family)) family=payload_family
    do p=1,nnode
       d=node_domain(p)
       id=node_id(p)
       do k=1,zlevels
          do v=first,last
             if (component==1) then
                var=scalars(1)+v-1
                address=id+1
             else
                var=S_VELO
                address=EDGE*id+v-nscalar
             end if
             native=value(v,k,p,family)
             reference=scaling(var,k)%data(d)%elts(address)
             if (transfer_bits(native)==transfer_bits(reference)) cycle
             if (.not. ieee_is_finite(native) .or. .not. ieee_is_finite(reference)) call die('nonfinite oracle difference')
             allowed=128.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(native),abs(reference))
             if (abs(native-reference) <= allowed) cycle
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

  subroutine publish_native_inverse(scaling,payload_family)
    ! One compatibility publication after the complete native transaction.
    ! Interior installation remains the standard final-owner publication.
    type(Float_Field), intent(inout) :: scaling(1:N_VARIABLE,1:zlevels)
    integer, optional, intent(in) :: payload_family
    integer :: p,d,id,k,v,family
    family=BLOCK_PAYLOAD_SOL
    if (present(payload_family)) family=payload_family
    do p=1,nnode
       d=node_domain(p)
       id=node_id(p)
       if (nodes(d)%key(1,id+1)>0) cycle
       do k=1,zlevels
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
             coverage(1,l)=coverage(1,l)+int(zlevels,int64)
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
                      coverage(2,l)=coverage(2,l)+int(zlevels,int64)
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
                         coverage(3,l)=coverage(3,l)+2_int64*zlevels
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
      if (item%kind/=0) coverage(4,l)=coverage(4,l)+2_int64*zlevels
    end subroutine
  end subroutine

  subroutine native_inverse_outer(level,count)
    integer, intent(in) :: level
    integer(int64), intent(out) :: count(4)
    integer :: p,k,q
    real(dp) :: s(13),x(9),correction,second,first
    count=coverage(:,level)
    first=0.0_dp
    second=0.0_dp
    do k=1,zlevels
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
  end subroutine
end module parallel_block_inverse_mod
