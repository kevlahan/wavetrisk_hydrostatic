module parallel_block_rhs_boundary_mod
  ! Topology-only translation of the accepted boundary transaction. Remote
  ! sends sample the pre-copy state; same-rank copies execute in original order;
  ! remote receives overwrite last. Resolve that ordering symbolically once.
  ! The resulting physical identities can be placed on arbitrary block owners.
  use kind_mod, only : dp
  use mpi_f08
  use arch_mod, only : rank,n_process,comm,glo_id
  use shared_mod, only : n_domain,AT_NODE,AT_EDGE,EDGE,NONE,level_start,level_end
  use domain_mod, only : grid
  use parallel_block_rhs_exchange_mod, only : Rhs_Exchange,compile_rhs_exchange
  use iso_fortran_env, only : int64
  implicit none
  private
  public :: Rhs_Boundary_Map,compile_rhs_boundary_map,bind_rhs_boundary_exchange
  type :: Rhs_Boundary_Map
     integer, allocatable :: key(:,:),source(:,:),target(:),orientation(:)
  end type
contains
  logical function matches_level(node_level,requested) result(matches)
    integer,intent(in) :: node_level,requested
    matches=node_level==requested
    if(requested==NONE) matches=node_level>=level_start-1.and.node_level<=level_end
  end function
  subroutine compile_rhs_boundary_map(pos,level,base,output)
    integer, intent(in) :: pos,level,base(:)
    type(Rhs_Boundary_Map), intent(out) :: output
    integer, allocatable :: sc(:),sd(:),rc(:),rd(:),sent(:,:),received(:,:),dest(:),sgn(:),map(:,:),signs(:)
    logical, allocatable :: changed(:)
    integer :: mult,n,ds,dd,r,g,i,id,ns,nr,pass,si,di,s,ierr
    if(pos/=AT_NODE.and.pos/=AT_EDGE) error stop 'RHS boundary location'
    if(size(base)/=size(grid)+1) error stop 'RHS boundary base extent'
    mult=1
    if(pos==AT_EDGE) mult=EDGE
    n=0
    do ds=1,size(grid)
       if(base(ds)/=n) error stop 'RHS boundary noncontiguous namespace'
       n=n+grid(ds)%node%length
    end do
    if(base(size(base))/=n) error stop 'RHS boundary total extent'
    n=mult*n
    allocate(output%key(4,n),map(4,n),signs(n),changed(n))
    do ds=1,size(grid)
       do id=0,mult*grid(ds)%node%length-1
          output%key(:,mult*base(ds)+id+1)=[glo_id(rank+1,ds)+1,id/mult+1,mod(id,mult)+1,1]
       end do
    end do
    map=output%key;signs=1;changed=.false.
    allocate(sc(n_process),sd(n_process),rc(n_process),rd(n_process))
    sc=0;rc=0
    do pass=1,2
       ns=0;nr=0
       do r=1,n_process
          sd(r)=ns;rd(r)=nr
          if(r==rank+1) cycle
          do ds=1,size(grid)
             do dd=1,n_domain(r)
                g=glo_id(r,dd)+1
                do i=1,grid(ds)%pack(pos,g)%length
                   id=grid(ds)%pack(pos,g)%elts(i)
                   if(.not.matches_level(grid(ds)%level%elts(id/mult+1),level)) cycle
                   ns=ns+1
                   if(pass==2) sent(:,ns)=output%key(:,mult*base(ds)+id+1)
                end do
             end do
          end do
          sc(r)=ns-sd(r)
          do ds=1,n_domain(r)
             g=glo_id(r,ds)+1
             do dd=1,size(grid)
                do i=1,grid(dd)%unpk(pos,g)%length
                   id=grid(dd)%unpk(pos,g)%elts(i)
                   if(.not.matches_level(grid(dd)%level%elts(abs(id)/mult+1),level)) cycle
                   nr=nr+1
                   if(pass==2) then
                      dest(nr)=mult*base(dd)+abs(id)+1
                      sgn(nr)=merge(-1,1,id<0.and.pos==AT_EDGE)
                   end if
                end do
             end do
          end do
          rc(r)=nr-rd(r)
       end do
       if(pass==1) allocate(sent(4,ns),received(4,nr),dest(nr),sgn(nr))
    end do
    call MPI_Alltoallv(sent,4*sc,4*sd,MPI_INTEGER,received,4*rc,4*rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'RHS boundary metadata delivery'
    ! cp_bdry_inside intentionally visits ALL levels, even in a level halo.
    ! Propagate source identity/sign after each write; a snapshot of every
    ! local source would silently change chains and cycles between aliases.
    do ds=1,size(grid)
       do dd=1,size(grid)
          g=glo_id(rank+1,dd)+1
          do i=1,grid(ds)%pack(pos,g)%length
             si=mult*base(ds)+grid(ds)%pack(pos,g)%elts(i)+1
             id=grid(dd)%unpk(pos,glo_id(rank+1,ds)+1)%elts(i)
             di=mult*base(dd)+abs(id)+1
             map(:,di)=map(:,si);signs(di)=merge(-1,1,id<0.and.pos==AT_EDGE)*signs(si)
             changed(di)=.true.
          end do
       end do
    end do
    do i=1,nr
       map(:,dest(i))=received(:,i);signs(dest(i))=sgn(i);changed(dest(i))=.true.
    end do
    allocate(output%target(count(changed)),output%source(4,count(changed)),output%orientation(count(changed)))
    s=0
    do i=1,n
       if(.not.changed(i)) cycle
       s=s+1;output%target(s)=i;output%source(:,s)=map(:,i);output%orientation(s)=signs(i)
    end do
  end subroutine

  subroutine bind_rhs_boundary_exchange(map,plan,topology,masks,phase,nlayer)
    ! Same-placement binding for qualification. Distributed placement uses
    ! these exact source/target identities with the block slot directory.
    type(Rhs_Boundary_Map), intent(in) :: map
    type(Rhs_Exchange), intent(inout) :: plan
    integer(int64), intent(in) :: topology,masks
    integer, intent(in) :: phase,nlayer
    integer, allocatable :: order(:,:),slot(:)
    integer :: i,n
    n=size(map%key,2);allocate(order(2,n),slot(n));order=1
    do i=1,n
       slot(i)=i
    end do
    call compile_rhs_exchange(plan,comm,topology,masks,phase,n,nlayer,map%key,order,slot, &
         map%source,map%target,map%orientation)
  end subroutine
end module
