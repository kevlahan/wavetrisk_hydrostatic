program test_inverse_routes
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : EDGE, Coord, zlevels
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod
  implicit none
  integer :: b,p,c,f,n,nv,first,e
  integer, allocatable :: seen(:),address(:,:),empty(:,:)
  integer :: keys(4,6)
  real(dp), allocatable :: wire(:),reference(:,:),roundtrip(:)

  real(dp) :: workspace(5,3,6)
  integer :: self_slots(2)=[5,2],selected(6)=[1,2,3,1,3,2]
  zlevels=3
  allocate(block_source(2),block_received(0),block_received_catalog_index(0))
  block_source_catalog_index=[2,1]
  block_retained_source_index=[1,2]
  do b=1,2
     call initialize(block_source(b),b)
  end do
  call install_local_blocks(2,seen)
  ! Mixed catalogs, disjoint storage, nonzero patch/boundary offsets and
  ! atmospheric layers embedded in a wider field-level extent.
  keys(:,1)=[2,STORE_PATCH,2,PATCH_SIZE**2-1]
  keys(:,2)=[1,STORE_BDRY,2,1]
  keys(:,3)=[2,STORE_BDRY,1,0]
  keys(:,4)=[1,STORE_PATCH,1,0]
  keys(:,5)=[2,STORE_PATCH,1,2]
  keys(:,6)=[1,STORE_BDRY,1,2]
  call compile_local_block_inverse_routes(keys,address)
  do f=BLOCK_PAYLOAD_SOL,BLOCK_PAYLOAD_WAV_COEFF
     do c=1,2
        nv=merge(2,EDGE,c==1)
        n=nv*zlevels
        allocate(wire(6*n),reference(nv,zlevels),roundtrip(6*n))
        call transfer_local_block_inverse_routes(address,f,c,.false.,wire)
        do p=1,6
           call transfer_local_block_inverse_node(keys(:,p),f,c,.false.,reference)
           call assert_equal(wire((p-1)*n+1:p*n),reshape(reference,[n]))
        end do
        wire=-wire-37.0_dp
        call transfer_local_block_inverse_routes(address,f,c,.true.,wire)
        do p=1,6
           call transfer_local_block_inverse_node(keys(:,p),f,c,.false.,reference)
           call assert_equal(wire((p-1)*n+1:p*n),reshape(reference,[n]))
           reference=reference+13.0_dp
           call transfer_local_block_inverse_node(keys(:,p),f,c,.true.,reference)
        end do
        call transfer_local_block_inverse_routes(address,f,c,.false.,roundtrip)
        call assert_equal(roundtrip,wire+13.0_dp)
        ! Treat the middle two routes as locally owned. Their wire slots must
        ! stay untouched; data travels straight between workspace and fields.
        first=merge(1,3,c==1)
        workspace=-999.0_dp
        wire=-777.0_dp
        call transfer_local_block_inverse_routes(address,f,c,.false.,wire, &
             workspace=workspace,local_first=3,local_slots=self_slots)
        call assert_equal(wire(1:2*n),roundtrip(1:2*n))
        call assert_equal(wire(4*n+1:),roundtrip(4*n+1:))
        if(any(wire(2*n+1:4*n)/=-777.0_dp))error stop 'local gather staged on wire'
        do p=1,2
           call assert_equal(reshape(workspace(first:first+nv-1,:,self_slots(p)),[n]),roundtrip((p+1)*n+1:(p+2)*n))
        end do
        workspace=workspace+7.0_dp
        wire=wire+7.0_dp
        call transfer_local_block_inverse_routes(address,f,c,.true.,wire, &
             workspace=workspace,local_first=3,local_slots=self_slots)
        call transfer_local_block_inverse_routes(address,f,c,.false.,wire)
        call assert_equal(wire,roundtrip+7.0_dp)
        deallocate(wire,reference,roundtrip)
     end do
  end do
  ! Selected-edge self installs preserve untouched vector components.
  do f=BLOCK_PAYLOAD_SOL,BLOCK_PAYLOAD_WAV_COEFF
     allocate(wire(6*zlevels),roundtrip(6*EDGE*zlevels),reference(EDGE,zlevels))
     call transfer_local_block_inverse_routes(address,f,2,.false.,roundtrip)
     workspace=-0.0_dp
     do p=1,6
        wire((p-1)*zlevels+1:p*zlevels)=real(-100-p,dp)
     end do
     workspace(5,:,5)=-303.0_dp
     workspace(3,:,2)=-404.0_dp
     call transfer_local_block_inverse_routes(address,f,2,.true.,wire,selected,workspace,3,self_slots)
     do p=1,6
        call transfer_local_block_inverse_node(keys(:,p),f,2,.false.,reference)
        do e=1,EDGE
           if(e==selected(p))then
              if(p==3.or.p==4)then
                 call assert_equal(reference(e,:),workspace(2+e,:,self_slots(p-2)))
              else
                 call assert_equal(reference(e,:),wire((p-1)*zlevels+1:p*zlevels))
              endif
           else
              call assert_equal(reference(e,:),roundtrip((p-1)*EDGE*zlevels+e:p*EDGE*zlevels:EDGE))
           endif
        enddo
     enddo
     deallocate(wire,roundtrip,reference)
  enddo
  allocate(empty(4,0),wire(0))
  call compile_local_block_inverse_routes(empty,address)
  call transfer_local_block_inverse_routes(address,BLOCK_PAYLOAD_SOL,1,.false.,wire)
  call transfer_local_block_inverse_routes(address,BLOCK_PAYLOAD_SOL,2,.true.,wire)
  call transfer_local_block_inverse_routes(address,BLOCK_PAYLOAD_SOL,1,.false.,wire, &
       workspace=workspace,local_first=1,local_slots=self_slots(:0))
  print *, 'PASS: inverse route pack/install, both families/components, field offsets and empty plans'
contains
  subroutine assert_equal(a,b)
    real(dp), intent(in) :: a(:),b(:)
    if (any(transfer(a,[0_int64],size(a))/=transfer(b,[0_int64],size(b)))) error stop 'inverse route mismatch'
  end subroutine
  subroutine initialize(block,identity)
    type(Block_Data), intent(out) :: block
    integer, intent(in) :: identity
    integer :: i,np,ns,nv
    block%id=identity
    block%field_level=-identity
    block%n_field_level=7
    block%n_scalar_variable=2
    np=2*PATCH_SIZE**2
    ns=2*7*np
    nv=EDGE*7*np
    allocate(block%patch(2),block%node(np),block%bdry_node(6),block%bdry_storage(2))
    do i=1,2
       block%patch(i)%elts_start=(i-1)*PATCH_SIZE**2
       block%patch(i)%level=5
       block%patch(i)%children=0
       block%patch(i)%neigh=0
       block%patch(i)%active=1
       block%patch(i)%deleted=.false.
       block%bdry_storage(i)%local_start=3*(i-1)
       block%bdry_storage(i)%n_node=3
    end do
    block%node=Coord(0.0_dp,0.0_dp,0.0_dp)
    block%bdry_node=Coord(0.0_dp,0.0_dp,0.0_dp)
    block%scalar=[(real(i+10000*identity,dp),i=1,ns)]
    block%vector=[(real(i+20000*identity,dp),i=1,nv)]
    block%wavelet_scalar=-block%scalar
    block%wavelet_vector=-block%vector
    block%bdry_scalar=[(real(i+30000*identity,dp),i=1,2*7*6)]
    block%bdry_vector=[(real(i+40000*identity,dp),i=1,EDGE*7*6)]
    block%bdry_wavelet_scalar=-block%bdry_scalar
    block%bdry_wavelet_vector=-block%bdry_vector
    allocate(block%scalar_mean(0),block%vector_mean(0),block%tke(0),block%wavelet_tke(0),block%topography(0))
    allocate(block%bdry_scalar_mean(0),block%bdry_vector_mean(0),block%bdry_tke(0), &
         block%bdry_wavelet_tke(0),block%bdry_topography(0))
    allocate(block%ghost_storage(0),block%ghost_node(0),block%ghost_scalar(0),block%ghost_vector(0), &
         block%ghost_wavelet_scalar(0),block%ghost_wavelet_vector(0),block%ghost_scalar_mean(0), &
         block%ghost_vector_mean(0),block%ghost_tke(0),block%ghost_wavelet_tke(0),block%ghost_topography(0))
    allocate(block%neigh_class(0,0),block%block_bdry(0),block%stencil(0,0))
  end subroutine
end program
