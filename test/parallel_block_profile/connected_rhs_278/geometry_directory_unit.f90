program geometry_directory_unit
  use iso_fortran_env,only:int64
  use kind_mod,only:dp
  use mpi_f08
  use arch_mod,only:comm,rank,n_process,glo_id
  use shared_mod,only:N_GLO_DOMAIN,n_domain
  use domain_mod,only:grid
  use parallel_block_rhs_distribute_mod
  implicit none
  type(Rhs_Geometry_Directory) :: directory
  integer :: ierr,epoch,owners,g,r,d,j,c,m,n,nq,i,k,width
  integer(int64) :: topology
  integer,allocatable :: base(:),key(:,:),query(:,:)
  real(dp),allocatable :: source(:,:),actual(:,:),expected(:,:)
  character(32) :: mode
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,n_process,ierr)
  call get_command_argument(1,mode)
  owners=max(1,n_process-1)
  allocate(n_domain(n_process))
  do epoch=1,12
     ! Change owner, local Domain order and node extents, including a rank
     ! without geometry which still requests metadata on every iteration.
     topology=int((epoch+1)/2,int64)
     n_domain=0
     do g=N_GLO_DOMAIN,1,-1
        r=1+modulo(g+epoch,owners);n_domain(r)=n_domain(r)+1
     end do
     allocate(glo_id(n_process,maxval(n_domain)),grid(n_domain(rank+1)),base(n_domain(rank+1)+1))
     glo_id=-1;n_domain=0
     do g=N_GLO_DOMAIN,1,-1
        r=1+modulo(g+epoch,owners);n_domain(r)=n_domain(r)+1;glo_id(r,n_domain(r))=g-1
     end do
     base(1)=0
     do d=1,size(grid)
        g=glo_id(rank+1,d)+1;grid(d)%node%length=2+modulo(g+epoch,5)
        base(d+1)=base(d)+grid(d)%node%length
     end do
     call prepare_rhs_geometry_directory(base,topology,directory)
     do m=1,3,2
        n=m*base(size(base));nq=13;width=1+modulo(epoch,4)
        if(modulo(epoch,2)==0.and.rank==0) nq=0
        allocate(key(4,n),source(n,width),query(4,nq))
        do d=1,size(grid)
           g=glo_id(rank+1,d)+1
           do j=1,grid(d)%node%length
              do c=1,m
                 i=m*(base(d)+j-1)+c;key(:,i)=[g,j,c,1]
                 do k=1,width
                    source(i,k)=real(100000*epoch+1000*g+100*j+10*c+k,dp)
                 end do
              end do
           end do
        end do
        do i=1,nq
           g=1+modulo(i/2+rank+epoch,N_GLO_DOMAIN)
           query(:,i)=[g,1+modulo(i,2+modulo(g+epoch,5)),1+modulo(i,m),1]
        end do
        if(trim(mode)=='source-key'.and.rank==0) key(4,1)=2
        if(trim(mode)=='request-key'.and.rank==0) query(4,1)=2
        if(trim(mode)=='layout'.and.rank==0) grid(1)%node%length=grid(1)%node%length+1
        if(trim(mode)=='stale') then
           call deliver_rhs_metadata(key,source,query,actual,topology+1,100+epoch,directory,m)
        else
           call deliver_rhs_metadata(key,source,query,actual,topology,100+epoch,directory,m)
        end if
        call deliver_rhs_metadata(key,source,query,expected,topology,100+epoch)
        if(any(transfer(actual,[0_int64],nq*width)/=transfer(expected,[0_int64],nq*width))) &
             error stop 'geometry directory differs from generic compiler'
        deallocate(key,source,query,actual,expected)
     end do
     deallocate(glo_id,grid,base)
  end do
  if(rank==0) print *, 'PASS geometry directory, ownership migration, node/component identities and invalidation'
  call MPI_Finalize(ierr)
end program
