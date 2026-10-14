program layout_unit
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : EDGE,N_BDRY,S_MASS,S_TEMP,S_VELO,EAST,WEST
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data
  use parallel_block_rhs_layout_mod
  implicit none
  type(Block_Data) :: b
  type(Rhs_Input_Layout) :: plan
  integer :: n,nb,p,s,i,k,j,classes(16)
  real(dp) :: x(16),before(16)
  logical :: rebuilt
  character(len=32) :: mode
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
  call compile_rhs_input_layout(b,1_int64,plan,rebuilt)
  if (.not.rebuilt) error stop 'initial layout not built'
  select case(trim(mode))
  case('stale')
     call rhs_input_values(b,plan,2_int64,1,1,1,x,classes)
  case('topology')
     b%patch(2)%neigh(EAST)=-EAST
     call compile_rhs_input_layout(b,1_int64,plan,rebuilt)
  case('missing')
     b%patch(1)%neigh(EAST)=0
     call compile_rhs_input_layout(b,2_int64,plan,rebuilt)
  case('layer')
     call rhs_input_values(b,plan,1_int64,1,1,3,x,classes)
  case('extent')
     deallocate(b%scalar);allocate(b%scalar(1))
     call compile_rhs_input_layout(b,2_int64,plan,rebuilt)
  case('')
     do k=0,2
        do p=1,3
           do j=0,PATCH_SIZE-1
              do i=0,PATCH_SIZE-1
                 call rhs_input_values(b,plan,1_int64,p,j*PATCH_SIZE+i+1,k,x,classes)
                 if (x(1)/=real(k*n+(p-1)*PATCH_SIZE**2+j*PATCH_SIZE+i+1,dp)+0.25_dp) &
                      error stop 'owned value mismatch'
                 if (classes(1)/=1) error stop 'owned storage mismatch'
              end do
           end do
        end do
     end do
     call rhs_input_values(b,plan,1_int64,2,PATCH_SIZE,1,before,classes)
     if (classes(3)/=1) error stop 'internal neighbor not owned'
     call rhs_input_values(b,plan,1_int64,1,PATCH_SIZE,1,x,classes)
     if (classes(3)/=3) error stop 'external neighbor not ghost'
     b%patch(2)%active=1
     call compile_rhs_input_layout(b,1_int64,plan,rebuilt)
     if (rebuilt) error stop 'mask-only change rebuilt input addresses'
     ! Masks are deliberately outside this input-layout contract; no numerical plan reuse is implied.
     b%patch(2)%neigh(EAST)=-EAST
     call compile_rhs_input_layout(b,2_int64,plan,rebuilt)
     if (.not.rebuilt) error stop 'new topology not rebuilt'
     call rhs_input_values(b,plan,2_int64,2,PATCH_SIZE,1,x,classes)
     if (classes(3)/=3.or.x(3)==before(3)) error stop 'changed topology not reflected'
     if (rhs_input_layout_bytes(plan)<=0_int64) error stop 'missing storage account'
     print *, 'PASS layout owned/ghost/boundary extents, reuse, invalidation and layers'
  case default
     error stop 'unknown test'
  end select
  if (len_trim(mode)>0) error stop 'guard did not reject'
end program
