program boundary_contract215
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : n_domain,AT_NODE,AT_EDGE,EDGE
  use arch_mod, only : init_arch_mod,comm,rank,n_process,glo_id,n_glo_block
  use domain_mod, only : grid,Float_Field
  use comm_mpi_mod, only : init_comm_mpi_mod,update_bdry1,req, &
       recv_lengths,recv_offsets,send_lengths,send_offsets
  use mpi_f08
  implicit none
  type(Float_Field) :: fields(3,2)
  integer :: d,g,p,k,i,j,r,mode,window,lo,hi,unit,counter
  character(100) :: name,arg
  call init_arch_mod
  call get_command_argument(1,arg)
  n_domain=2
  ! Exercise a genuinely idle communication participant, without claiming
  ! that the full climate startup supports more ranks than geometric Domains.
  if(n_process>1) n_domain(n_process)=0
  n_glo_block=sum(n_domain)
  allocate(glo_id(n_process,2));glo_id=-1;g=0
  do r=1,n_process
     do d=1,n_domain(r)
        glo_id(r,d)=g;g=g+1
     end do
  end do
  allocate(grid(n_domain(rank+1)))
  allocate(req(2*n_process),recv_lengths(n_process),recv_offsets(n_process), &
       send_lengths(n_process),send_offsets(n_process))
  call init_comm_mpi_mod
  do d=1,size(grid)
     allocate(grid(d)%level%elts(32));grid(d)%level%length=32
     grid(d)%level%elts=[(mod(i-1,3)+1,i=1,32)]
     allocate(grid(d)%pack(AT_NODE:AT_EDGE,n_glo_block),grid(d)%unpk(AT_NODE:AT_EDGE,n_glo_block))
     do g=1,n_glo_block
        do p=AT_NODE,AT_EDGE
           allocate(grid(d)%pack(p,g)%elts(4),grid(d)%unpk(p,g)%elts(4))
           grid(d)%pack(p,g)%length=4;grid(d)%unpk(p,g)%length=4
           if(p==AT_NODE) then
              grid(d)%pack(p,g)%elts=[0,1,2,3]
              grid(d)%unpk(p,g)%elts=[3,4,5,0]
           else
              grid(d)%pack(p,g)%elts=[0,3,6,9]
              grid(d)%unpk(p,g)%elts=[-9,12,-15,0]
           end if
        end do
     end do
  end do
  do j=1,2
     do i=1,3
        fields(i,j)%pos=merge(AT_EDGE,AT_NODE,i==2)
        allocate(fields(i,j)%data(size(grid)))
        do d=1,size(grid)
           allocate(fields(i,j)%data(d)%elts(96));fields(i,j)%data(d)%length=96
        end do
     end do
  end do
  write(name,'(a,i0,a)') 'boundary-',rank,'.txt'
  open(newunit=unit,file=trim(name),status='replace');counter=0
  do mode=0,2
     do window=1,5
        lo=1;hi=3
        if(window==2.or.window==4) lo=2
        if(window==2.or.window==4) hi=2
        if(window==5) lo=4
        if(window==5) hi=4
        do j=1,2
           do i=1,3
              do d=1,size(grid)
                 fields(i,j)%data(d)%elts=[(real(100000*rank+10000*j+1000*i+100*d+k,dp),k=1,96)]
                 fields(i,j)%data(d)%elts(1)=0.0_dp
              end do
           end do
        end do
        fields%bdry_uptodate=.false.;fields%bdry_tag=-1
        if(mode==0) then
           do j=1,2
              do i=1,3
                 call update_bdry1(fields(i,j),lo,hi)
              end do
           end do
        else if(mode==1) then
           do j=1,2
              call update_bdry1(fields(:,j),lo,hi)
           end do
        else
           call update_bdry1(fields,lo,hi)
        end if
        ! An unchanged request with all fields current must be a no-op.
        if(lo<hi) call update_bdry1(fields,lo,hi)
        do j=1,2
           do i=1,3
              do d=1,size(grid)
                 do k=1,96
                    write(unit,'(z16.16)') transfer(fields(i,j)%data(d)%elts(k),0_int64)
                    counter=counter+1
                 end do
              end do
           end do
        end do
     end do
  end do
  close(unit)
  if(trim(arg)=='stale') then
     ! Keep list lengths AND level counts unchanged; only the address changes.
     ! The checked candidate must reject it before posting the next exchange.
     do d=1,size(grid)
        do g=1,n_glo_block
           grid(d)%pack(AT_NODE,g)%elts(1)=3
        end do
     end do
     fields%bdry_uptodate=.false.
     call update_bdry1(fields,1,3)
  end if
  call MPI_Finalize
end program
