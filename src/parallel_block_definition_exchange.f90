module parallel_block_definition_exchange_mod
  ! Version/address traffic only. Preserve update_bdry's remote snapshot,
  ! ordered all-level local copies, then remote installation.
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,glo_id,abort_run
  use shared_mod, only: n_domain,AT_EDGE
  use domain_mod, only: grid
  implicit none
  private
  public :: Owner_Definitions,exchange_edge_definitions
  type :: Owner_Definitions
     integer,allocatable :: key(:,:)
  end type
contains
  subroutine exchange_edge_definitions(field,l,coverage)
    type(Owner_Definitions),intent(inout)::field(:)
    integer,intent(in)::l
    logical,intent(in),optional::coverage
    logical::mark
    integer::sc(n_process),sd(n_process),rc(n_process),rd(n_process), &
         r,ds,dd,g,i,id,n,a,b,ns,nr,pass,ierr,orientation
    integer,allocatable::send(:,:),recv(:,:),target(:,:)
    integer::value(7)
    mark=.false.
    if(present(coverage)) mark=coverage
    do pass=1,2
       ns=0;nr=0;sc=0;rc=0
       do r=1,n_process
          sd(r)=ns;rd(r)=nr
          if(r==rank+1) cycle
          do ds=1,size(grid)
             do dd=1,n_domain(r)
                g=glo_id(r,dd)+1
                do i=1,grid(ds)%pack(AT_EDGE,g)%length
                   id=grid(ds)%pack(AT_EDGE,g)%elts(i)
                   if(grid(ds)%level%elts(id/3+1)/=l) cycle
                   ns=ns+1
                   if(pass==2) send(:,ns)=field(ds)%key(:,id+1)
                end do
             end do
          end do
          do ds=1,n_domain(r)
             g=glo_id(r,ds)+1
             do dd=1,size(grid)
                do i=1,grid(dd)%unpk(AT_EDGE,g)%length
                   id=grid(dd)%unpk(AT_EDGE,g)%elts(i)
                   if(grid(dd)%level%elts(abs(id)/3+1)/=l) cycle
                   nr=nr+1
                   if(pass==2) target(:,nr)=[dd,abs(id)+1,merge(-1,1,id<0)]
                end do
             end do
          end do
          sc(r)=ns-sd(r);rc(r)=nr-rd(r)
       end do
       if(pass==1) allocate(send(7,ns),recv(7,nr),target(3,nr))
    end do
    call MPI_Alltoallv(send,7*sc,7*sd,MPI_INTEGER,recv,7*rc,7*rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call abort_run
    do ds=1,size(grid)
       do dd=1,size(grid)
          g=glo_id(rank+1,dd)+1
          do i=1,grid(ds)%pack(AT_EDGE,g)%length
             a=grid(ds)%pack(AT_EDGE,g)%elts(i)+1
             id=grid(dd)%unpk(AT_EDGE,glo_id(rank+1,ds)+1)%elts(i);b=abs(id)+1
             orientation=merge(-1,1,id<0)
             value=field(ds)%key(:,a);value(5)=orientation*value(5)
             if(mark) value(1)=99
             field(dd)%key(:,b)=value
          end do
       end do
    end do
    do n=1,nr
       dd=target(1,n);b=target(2,n)
       value=recv(:,n);value(5)=target(3,n)*value(5)
       if(mark) value(1)=99
       field(dd)%key(:,b)=value
    end do
  end subroutine
end module
