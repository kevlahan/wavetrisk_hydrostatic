program schedule_unit
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : EDGE,N_BDRY,S_MASS,S_TEMP,S_VELO,EAST,WEST,SOUTH,SOUTHEAST,IJPLUS,IJMINUS,IPLUSJMINUS,IMINUSJPLUS
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data
  use rhs_schedule_check_mod
  use parallel_block_rhs_schedule_mod
  use parallel_block_rhs_layout_mod
  implicit none
  type(Block_Data) :: b
  type(Rhs_Input_Layout) :: plan
  integer :: n,nb,p,s,i,k,j,classes(16)
  real(dp) :: x(16),before(16)
  logical :: rebuilt
  character(len=32) :: mode
  type(Rhs_Producer_Schedule) :: program
  integer(int64) :: counts(5)
  integer :: flags,corners(4)
  integer, allocatable :: fw(:,:),pvw(:,:),vw(:,:),keys(:,:),dk(:,:),orders(:,:),slots(:)
  integer :: f,m,w,q,ordinal(3),within,seen
  call get_command_argument(1,mode)
  n=3*PATCH_SIZE**2;nb=N_BDRY*PATCH_SIZE**2
  b%id=1;b%root_domain=0;b%root_patch=0;b%scalar_variable=S_MASS
  b%n_scalar_variable=S_TEMP-S_MASS+1;b%vector_variable=S_VELO
  b%field_level=0;b%n_field_level=3;b%scalar_mult=1;b%vector_mult=EDGE
  allocate(b%patch(3),b%node(n),b%bdry_node(nb),b%ghost_node(PATCH_SIZE**2))
  allocate(b%scalar(6*n),b%scalar_mean(6*n),b%vector(9*n))
  allocate(b%bdry_scalar(6*nb),b%bdry_scalar_mean(6*nb),b%bdry_vector(9*nb))
  allocate(b%ghost_scalar(6*PATCH_SIZE**2),b%ghost_scalar_mean(6*PATCH_SIZE**2),b%ghost_vector(9*PATCH_SIZE**2))
  allocate(b%block_bdry(N_BDRY),b%bdry_storage(N_BDRY),b%ghost_storage(1))
  do p=1,3
     b%patch(p)%elts_start=(p-1)*PATCH_SIZE**2;b%patch(p)%level=2
     b%patch(p)%deleted=.false.;b%patch(p)%active=0;b%patch(p)%children=0
     do s=1,N_BDRY
        b%patch(p)%neigh(s)=-s
     end do
  end do
  b%patch(2)%neigh(EAST)=2;b%patch(3)%neigh(WEST)=1
  do s=1,N_BDRY
     b%block_bdry(s)%storage_id=s;b%block_bdry(s)%ghost_id=0
     b%block_bdry(s)%dims=PATCH_SIZE;b%bdry_storage(s)%local_start=(s-1)*PATCH_SIZE**2
  end do
  b%block_bdry(EAST)%ghost_id=1;b%ghost_storage(1)%local_start=0
  b%ghost_storage(1)%n_node=PATCH_SIZE**2
  b%ghost_storage(1)%source_domain=1;b%ghost_storage(1)%source_block=2
  b%ghost_storage(1)%source_block_id=2;b%ghost_storage(1)%source_owner=1
  b%ghost_storage(1)%source_local_patch=1;b%ghost_storage(1)%source_level=2
  do i=1,size(b%scalar)
     b%scalar(i)=real(i,dp);b%scalar_mean(i)=0.25_dp
  end do
  do i=1,size(b%vector)
     b%vector(i)=real(10000+i,dp)
  end do
  do i=1,size(b%bdry_scalar)
     b%bdry_scalar(i)=real(20000+i,dp);b%bdry_scalar_mean(i)=0.5_dp
  end do
  do i=1,size(b%bdry_vector)
     b%bdry_vector(i)=real(30000+i,dp)
  end do
  do i=1,size(b%ghost_scalar)
     b%ghost_scalar(i)=real(40000+i,dp);b%ghost_scalar_mean(i)=0.75_dp
  end do
  do i=1,size(b%ghost_vector)
     b%ghost_vector(i)=real(50000+i,dp)
  end do
  corners=[IJPLUS,IPLUSJMINUS,IJMINUS,IMINUSJPLUS]
  call check_schedule(b,1_int64,counts)
  b%block_bdry(EAST)%ghost_id=0
  b%block_bdry(SOUTH)%dims=[PATCH_SIZE,3]
  b%block_bdry(SOUTHEAST)%dims=[3,PATCH_SIZE]
  b%block_bdry(EAST)%dims=[3,PATCH_SIZE]
  b%block_bdry%bdry_side=1
  call compile_rhs_input_layout(b,1_int64,plan,rebuilt)
  select case(trim(mode))
  case('stale')
     call compile_rhs_producer_schedule(b,plan,2_int64,2,[1,2,3],[integer::],[integer::],program)
  case('duplicate')
     call compile_rhs_producer_schedule(b,plan,1_int64,2,[1,1,2],[integer::],[integer::],program)
  case('incomplete')
     call compile_rhs_producer_schedule(b,plan,1_int64,2,[1,2],[integer::],[integer::],program)
  case('level')
     call compile_rhs_producer_schedule(b,plan,1_int64,3,[1,2,3],[integer::],[integer::],program)
  case('corner-order')
     call compile_rhs_producer_schedule(b,plan,1_int64,2,[1,2,3],[1,1],[IJMINUS,IJPLUS],program)
  case('corner-invalid')
     call compile_rhs_producer_schedule(b,plan,1_int64,2,[1,2,3],[1],[0],program)
  case('boundary-epoch')
     b%block_bdry(1)%bdry_side=0
     call compile_rhs_input_layout(b,1_int64,plan,rebuilt)
  case('missing-storage')
     b%block_bdry(1)%storage_id=0
     call compile_rhs_input_layout(b,2_int64,plan,rebuilt)
  case('')
     do flags=0,15
        do s=1,4
           b%block_bdry(s)%bdry_side=merge(1,0,btest(flags,s-1))
        end do
        call check_schedule(b,int(flags+2,int64),counts)
     end do
     b%block_bdry%bdry_side=1
     call compile_rhs_input_layout(b,18_int64,plan,rebuilt)
     call compile_rhs_producer_schedule(b,plan,18_int64,2,[3,2,1],[1,1,1,1],corners,program)
     allocate(fw(3,program%nslot),pvw(3,program%nslot),vw(2,program%nslot))
     call rhs_schedule_last_writes(program,fw,pvw,vw)
     ! All four mapped pentagon corrections own their final triangle writes.
     do s=1,4
        n=size(program%patch)+s
        select case(corners(s))
        case(IJPLUS)
           if (any(vw(:,program%pentagon_node(1,s))/=n)) error stop 'pentagon NE writer'
           if (pvw(1,program%pentagon_node(1,s))/=n.or.pvw(3,program%pentagon_node(1,s))/=n) &
                error stop 'pentagon NE PV writer'
        case(IPLUSJMINUS)
           if (vw(1,program%pentagon_node(7,s))/=n.or.vw(2,program%pentagon_node(4,s))/=n) &
                error stop 'pentagon SE writer'
           if (pvw(1,program%pentagon_node(1,s))/=n.or.pvw(2,program%pentagon_node(7,s))/=n) &
                error stop 'pentagon SE PV writer'
        case(IJMINUS)
           if (any(vw(:,program%pentagon_node(7,s))/=n)) error stop 'pentagon SW writer'
           if (pvw(1,program%pentagon_node(5,s))/=n.or.pvw(3,program%pentagon_node(4,s))/=n) &
                error stop 'pentagon SW PV writer'
        case(IMINUSJPLUS)
           if (vw(2,program%pentagon_node(7,s))/=n.or.vw(1,program%pentagon_node(5,s))/=n) &
                error stop 'pentagon NW writer'
           if (pvw(3,program%pentagon_node(1,s))/=n.or.pvw(2,program%pentagon_node(7,s))/=n) &
                error stop 'pentagon NW PV writer'
        end select
     end do
     allocate(keys(3,program%nslot))
     do i=1,program%nslot
        keys(:,i)=[2,1,i]
     end do
     ! This block has noncontiguous ordinals within a larger Domain traversal.
     ! All pentagons must follow patch 19 on a different owner.
     ordinal=[12,7,3]
     do f=1,3
        call rhs_schedule_definitions(program,ordinal,19,keys,f,11,dk,orders,slots)
        m=3
        if (f==3) m=2
        seen=0
        do i=1,size(slots)
           p=(slots(i)-12)/m+1;c_loop: block
           integer :: c
           c=mod(slots(i)-12,m)+1
           if (any(dk(:,i)/=[keys(:,p),c])) error stop 'writer key or address differs'
           select case(f)
           case(1)
              w=fw(c,p)
           case(2)
              w=pvw(c,p)
           case(3)
              w=vw(c,p)
           end select
           end block c_loop
           if (w==0) error stop 'undefined writer exported'
           if (abs(w)<=size(program%patch)) then
              q=abs(w)
              within=count(program%patch(:q)==program%patch(q))
              if (any(orders(:,i)/=[ordinal(program%patch(q)),2*within-merge(1,0,w>0)])) &
                   error stop 'direct or SW writer order differs'
           else
              if (any(orders(:,i)/=[20,program%corner(w-size(program%patch))])) &
                   error stop 'pentagon global writer order differs'
           end if
           seen=seen+1
        end do
        select case(f)
        case(1)
           if (seen/=count(fw/=0)) error stop 'flux writer coverage'
        case(2)
           if (seen/=count(pvw/=0)) error stop 'PV writer coverage'
        case(3)
           if (seen/=count(vw/=0)) error stop 'vorticity writer coverage'
        end select
     end do
     if (rhs_schedule_bytes(program)<=0_int64) error stop 'schedule storage not accounted'
     print *, 'PASS schedule all physical boundary combinations, both patch orders, pentagon writers and guards'
  case default
     error stop 'unknown test'
  end select
  if (len_trim(mode)>0) error stop 'guard did not reject'
end program
