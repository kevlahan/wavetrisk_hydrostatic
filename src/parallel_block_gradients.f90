module parallel_block_gradients_mod
  use parallel_block_velocity_owner_mod, only: velocity_owner_enabled,owner_velocity_source_key,owner_velocity_source_value
  ! Final velocity gradients consume the actual block primitive results after
  ! the compiled Bernoulli/Exner restrictions. The complete owner path reads
  ! source terms from their actual block producers; the earlier feature mode
  ! retains an explicit native source input. No q/B/Exner Domain input is read.
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run,block_catalog
  use shared_mod, only: zlevels
  use patch_mod, only: PATCH_SIZE
  use parallel_block_mod, only: local_block_patch_count,STORE_PATCH
  use parallel_block_inverse_mod, only: inverse_operand_generation,block_rk_enabled
  use parallel_block_primitives_mod, only: primitives_active,primitives_oracle, &
       primitive_thermodynamic_source,primitive_thermodynamic_value
  use parallel_block_resident_operands_mod
  use parallel_block_velocity_mod, only: native_velocity,velocity_gradient
  implicit none
  private
  public :: execute_block_gradients,install_block_gradients
  public :: owned_gradient_node
  public :: reset_owned_gradients,owned_gradients_ready,copy_owned_gradient_patch
  type :: Block_Gradient_Map
     integer, allocatable :: operation(:,:)
  end type
  type(Block_Gradient_Map), allocatable :: block_map(:)
  real(dp), allocatable :: retained(:,:,:)
  logical, allocatable :: layer_ready(:)
  logical :: owned_ready=.false.
  integer, allocatable :: sc(:),sd(:),rc(:),rd(:),identity(:,:),order(:)
  integer,allocatable :: vc(:),vd(:),wc(:),wd(:),velocity_request(:,:),velocity_address(:)
  real(dp),allocatable :: velocity_values(:)
  integer, allocatable :: tc(:),td(:),uc(:),ud(:),thermo_request(:,:),thermo_address(:)
  real(dp), allocatable :: length(:,:),source(:,:),computed(:,:),returned(:,:),thermo(:,:)
  type(Resident_Operand_Plan) :: input
  integer(int64) :: generation=-1_int64,executed=0_int64,checked=0_int64
  integer :: pending_layer=0
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': block gradients: ',message
    call abort_run
    error stop 'block gradient abort returned'
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
  subroutine compile_plan
    integer :: d,s,n,i,j,r,p,t,b,ierr,keys(5,4),producer(4),node(4),vkey(5)
    integer, allocatable :: count(:),send_meta(:,:),meta(:,:),dest(:),packed(:,:)
    real(dp), allocatable :: send_length(:,:)
    if(allocated(sc)) deallocate(sc,sd,rc,rd,identity,order,tc,td,uc,ud,thermo_request,thermo_address, &
         length,source,computed,returned,thermo,block_map,retained,layer_ready)
    if(allocated(vc)) deallocate(vc,vd,wc,wd,velocity_request,velocity_address,velocity_values)
    n=0
    do d=1,size(native_velocity)
       n=n+count_active(d)
    end do
    allocate(identity(2,n),order(n),dest(n),sc(n_process),sd(n_process),rc(n_process),rd(n_process),count(n_process))
    sc=0;i=0
    do d=1,size(native_velocity)
       do s=1,size(native_velocity(d)%gradient)
          if(.not.native_velocity(d)%gradient(s)%active) cycle
          i=i+1;identity(:,i)=[d,s]
          keys(:,1)=resident_operand_key(d,native_velocity(d)%gradient(s)%node(1)-1)
          if(keys(2,1)==0) call die('active gradient target has no compact block owner')
          dest(i)=keys(1,1);sc(dest(i)+1)=sc(dest(i)+1)+1
       end do
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('gradient counts')
    call offsets(sc,sd);call offsets(rc,rd)
    allocate(send_meta(51,n),meta(51,sum(rc)),send_length(3,n),length(3,sum(rc)))
    count=0
    do i=1,n
       d=identity(1,i);s=identity(2,i);node=native_velocity(d)%gradient(s)%node
       r=dest(i)+1;p=sd(r)+count(r)+1;count(r)=count(r)+1;order(i)=p
       do j=1,4
          keys(:,j)=resident_operand_key(d,node(j)-1)
          producer=primitive_thermodynamic_source(d,node(j))
          send_meta(20+4*j-3:20+4*j,p)=producer
       end do
       send_meta(37:51,p)=0
       if(velocity_owner_enabled()) then
          do j=1,3
             vkey=owner_velocity_source_key(d,node(1),j)
             send_meta(37+5*(j-1):36+5*j,p)=vkey
          end do
       end if
       send_meta(1:20,p)=reshape(keys,[20]);send_length(:,p)=native_velocity(d)%gradient(s)%length
    end do
    call MPI_Alltoallv(send_meta,51*sc,51*sd,MPI_INTEGER,meta,51*rc,51*rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('gradient input addresses')
    call MPI_Alltoallv(send_length,3*sc,3*sd,MPI_DOUBLE_PRECISION,length,3*rc,3*rd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('gradient geometry')
    allocate(packed(5,4*sum(rc)))
    do i=1,sum(rc)
       packed(:,4*i-3:4*i)=reshape(meta(1:20,i),[5,4])
    end do
    call compile_resident_operands(input,packed)
    deallocate(packed)
    allocate(tc(n_process),td(n_process),uc(n_process),ud(n_process),thermo_address(4*sum(rc)))
    tc=0
    do i=1,sum(rc)
       do j=1,4
          r=meta(20+4*j-3,i)+1
          if(r<1.or.r>n_process) call die('thermodynamic owner')
          tc(r)=tc(r)+1
       end do
    end do
    call MPI_Alltoall(tc,1,MPI_INTEGER,uc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('thermodynamic counts')
    call offsets(tc,td);call offsets(uc,ud)
    allocate(packed(3,sum(tc)),thermo_request(3,sum(uc)),thermo(2,sum(tc)))
    count=0;t=0
    do i=1,sum(rc)
       do j=1,4
          b=20+4*j-3;r=meta(b,i)+1;p=td(r)+count(r)+1;count(r)=count(r)+1
          packed(:,p)=meta(b+1:b+3,i);t=t+1;thermo_address(t)=p
       end do
    end do
    call MPI_Alltoallv(packed,3*tc,3*td,MPI_INTEGER,thermo_request,3*uc,3*ud,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('thermodynamic requests')
    if(velocity_owner_enabled()) then
       deallocate(packed)
       allocate(vc(n_process),vd(n_process),wc(n_process),wd(n_process),velocity_address(3*sum(rc)))
       vc=0
       do i=1,sum(rc)
          do j=1,3
             b=37+5*(j-1);r=meta(b,i)+1
             if(r<1.or.r>n_process) call die('source producer owner')
             vc(r)=vc(r)+1
          end do
       end do
       call MPI_Alltoall(vc,1,MPI_INTEGER,wc,1,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('source producer counts')
       call offsets(vc,vd);call offsets(wc,wd)
       allocate(packed(5,sum(vc)),velocity_request(5,sum(wc)),velocity_values(sum(vc)))
       count=0;t=0
       do i=1,sum(rc)
          do j=1,3
             b=37+5*(j-1);r=meta(b,i)+1;p=vd(r)+count(r)+1;count(r)=count(r)+1
             packed(:,p)=meta(b:b+4,i);t=t+1;velocity_address(t)=p
          end do
       end do
       call MPI_Alltoallv(packed,5*vc,5*vd,MPI_INTEGER,velocity_request,5*wc,5*wd,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('source producer requests')
    end if
    allocate(source(3,sum(rc)),computed(3,sum(rc)),returned(3,merge(0,n,block_rk_enabled.and..not.primitives_oracle)))
    allocate(block_map(size(block_catalog)),retained(3,sum(rc),zlevels),layer_ready(zlevels))
    layer_ready=.false.;owned_ready=.false.
    do b=1,size(block_catalog)
       if(block_catalog(b)%owner/=rank) cycle
       allocate(block_map(b)%operation(PATCH_SIZE**2,local_block_patch_count(b)))
       block_map(b)%operation=0
    end do
    do i=1,sum(rc)
       call register_retained_node(meta(2:5,i),i)
    end do
    call compile_inactive_targets
    do b=1,size(block_map)
       if(.not.allocated(block_map(b)%operation)) cycle
       if(any(block_map(b)%operation==0)) call die('retained block contains a node without a gradient definition')
    end do
    generation=inverse_operand_generation()
    write(*,'(a,i0,a,2(1x,i0))') 'Block gradient plan rank=',rank,' geometry/owner=',n,sum(rc)
  end subroutine
  subroutine register_retained_node(key,operation)
    integer, intent(in) :: key(4),operation
    integer :: b,p,node
    b=key(1);p=key(3);node=key(4)+1
    if(key(2)/=STORE_PATCH.or.b<1.or.b>size(block_map)) call die('retained gradient storage key')
    if(.not.allocated(block_map(b)%operation)) call die('retained gradient delivered to wrong owner')
    if(p<1.or.p>size(block_map(b)%operation,2).or.node<1.or.node>PATCH_SIZE**2) call die('retained node key')
    if(block_map(b)%operation(node,p)/=0) call die('duplicate retained gradient definition')
    block_map(b)%operation(node,p)=operation
  end subroutine
  subroutine compile_inactive_targets
    integer :: d,s,i,n,r,p,key(5),ierr,pass
    integer, allocatable :: counts(:),send_count(:),recv_count(:),send_disp(:),recv_disp(:),packed(:,:),received(:,:)
    allocate(counts(n_process),send_count(n_process),recv_count(n_process),send_disp(n_process),recv_disp(n_process))
    send_count=0
    do pass=1,2
       counts=0;n=0
       do d=1,size(native_velocity)
          do s=1,size(native_velocity(d)%gradient)
             if(native_velocity(d)%gradient(s)%active) cycle
             key=resident_operand_key(d,native_velocity(d)%gradient(s)%node(1)-1)
             if(key(2)==0) cycle ! non-block scaffold; only the RK return needs it
             r=key(1)+1;n=n+1
             if(pass==1) then
                send_count(r)=send_count(r)+1
             else
                p=send_disp(r)+counts(r)+1;counts(r)=counts(r)+1;packed(:,p)=key(2:5)
             end if
          end do
       end do
       if(pass==1) then
          call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,ierr)
          if(ierr/=MPI_SUCCESS) call die('inactive gradient counts')
          call offsets(send_count,send_disp);call offsets(recv_count,recv_disp)
          allocate(packed(4,n),received(4,sum(recv_count)))
       end if
    end do
    call MPI_Alltoallv(packed,4*send_count,4*send_disp,MPI_INTEGER, &
         received,4*recv_count,4*recv_disp,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('inactive gradient definitions')
    do i=1,size(received,2)
       call register_retained_node(received(:,i),-1)
    end do
  end subroutine
  subroutine reset_owned_gradients
    if(pending_layer/=0) call die('new transaction before gradient return consumed')
    owned_ready=.false.
    if(allocated(layer_ready)) layer_ready=.false.
  end subroutine
  logical function owned_gradients_ready()
    owned_gradients_ready=owned_ready
  end function
  function owned_gradient_node(b,p,node,k) result(x)
    integer,intent(in)::b,p,node,k
    real(dp)::x(3)
    integer::i
    if(.not.owned_ready.or.generation/=inverse_operand_generation()) call die('RK gradient epoch')
    if(.not.allocated(block_map(b)%operation)) call die('RK gradient non-owner')
    i=block_map(b)%operation(node+1,p)
    if(i==0) call die('RK gradient undefined')
    x=0.0_dp
    if(i>0) x=retained(:,i,k)
  end function

  subroutine copy_owned_gradient_patch(catalog,patch,first_level,nlevel,mult,values)
    integer, intent(in) :: catalog,patch,first_level,nlevel,mult ! patch is zero based
    real(dp), intent(out) :: values(:)
    integer :: j,k,slot,i,p
    if(.not.owned_ready.or.generation/=inverse_operand_generation()) call die('retained gradient epoch')
    if(mult/=3.or.size(values)/=3*PATCH_SIZE**2*nlevel) call die('retained gradient layout')
    if(catalog<1.or.catalog>size(block_map)) call die('retained gradient catalogue')
    if(.not.allocated(block_map(catalog)%operation)) call die('retained gradient owner')
    if(patch<0.or.patch>=size(block_map(catalog)%operation,2)) call die('retained gradient patch')
    values=0.0_dp
    do slot=1,nlevel
       k=first_level+slot-1
       if(k<1) cycle
       if(k>zlevels) call die('retained gradient layer')
       do j=1,PATCH_SIZE**2
          i=block_map(catalog)%operation(j,patch+1)
          if(i==0) call die('missing retained gradient definition')
          if(i<0) cycle ! explicit inactive definition
          p=3*((slot-1)*PATCH_SIZE**2+j-1)
          values(p+1:p+3)=retained(:,i,k)
       end do
    end do
  end subroutine
  integer function count_active(d) result(n)
    integer, intent(in) :: d
    integer :: s
    n=0
    do s=1,size(native_velocity(d)%gradient)
       if(native_velocity(d)%gradient(s)%active) n=n+1
    end do
  end function
  subroutine execute_block_gradients(k)
    integer, intent(in) :: k
    integer :: i,j,d,s,node,a,t,ierr
    real(dp), allocatable :: send_source(:,:),send_thermo(:,:),owner_source(:)
    real(dp) :: rho(4),theta_mass(4),b(4),e(4)
    if(.not.primitives_active) call die('gradient outside primitive transaction')
    if(pending_layer/=0) call die('unconsumed gradient results')
    if(generation/=inverse_operand_generation()) call compile_plan
    allocate(send_source(3,size(order)),send_thermo(2,sum(uc)))
    if(velocity_owner_enabled()) then
       allocate(owner_source(sum(wc)))
       do i=1,sum(wc)
          owner_source(i)=owner_velocity_source_value(velocity_request(:,i),k)
       end do
       call MPI_Alltoallv(owner_source,wc,wd,MPI_DOUBLE_PRECISION, &
            velocity_values,vc,vd,MPI_DOUBLE_PRECISION,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('owner source delivery')
       do i=1,size(source,2)
          do j=1,3
             source(j,i)=velocity_values(velocity_address(3*(i-1)+j))
          end do
       end do
    else
    do i=1,size(order)
       d=identity(1,i);s=identity(2,i);node=native_velocity(d)%gradient(s)%node(1)
       send_source(:,order(i))=native_velocity(d)%source(:,node)
    end do
    call MPI_Alltoallv(send_source,3*sc,3*sd,MPI_DOUBLE_PRECISION,source,3*rc,3*rd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('explicit source delivery')
    end if
    do i=1,sum(uc)
       send_thermo(:,i)=primitive_thermodynamic_value(thermo_request(:,i),k)
    end do
    call MPI_Alltoallv(send_thermo,2*uc,2*ud,MPI_DOUBLE_PRECISION,thermo,2*tc,2*td,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('thermodynamic result delivery')
    call fetch_resident_operands(input,k)
    do i=1,size(computed,2)
       do j=1,4
          t=4*(i-1)+j;a=input%address(t)
          rho(j)=input%recv(1,a);theta_mass(j)=input%recv(2,a)
          a=thermo_address(t);b(j)=thermo(1,a);e(j)=thermo(2,a)
       end do
       computed(:,i)=velocity_gradient(source(:,i),length(:,i),b,e,rho,theta_mass,.true.)
       executed=executed+1_int64
    end do
    ! Explicit return for the remaining RK closure, not a tendency operand.
    if(.not.block_rk_enabled.or.primitives_oracle) then
    call MPI_Alltoallv(computed,3*rc,3*rd,MPI_DOUBLE_PRECISION,returned,3*sc,3*sd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('gradient compatibility return')
    end if
    retained(:,:,k)=computed;layer_ready(k)=.true.;owned_ready=all(layer_ready)
    pending_layer=k
  end subroutine
  subroutine install_block_gradients(k)
    integer, intent(in) :: k
    integer :: i,d,s,node
    real(dp) :: actual(3),expected(3)
    if(pending_layer/=k) call die('gradient consumer epoch')
    if(block_rk_enabled.and..not.primitives_oracle) then
       pending_layer=0
       write(*,'(a,i0,a,2(1x,i0))') 'Block gradients rank=',rank,' executed/checked=',executed,checked
       return
    end if
    do i=1,size(order)
       d=identity(1,i);s=identity(2,i);node=native_velocity(d)%gradient(s)%node(1)
       expected=returned(:,order(i))
       if(primitives_oracle) then
          actual=native_velocity(d)%tendency(3*node-2:3*node,k)
          if(any(transfer(actual,[0_int64],3)/=transfer(expected,[0_int64],3))) then
             write(*,*) 'gradient mismatch d/node/layer',d,node,k,'actual',actual,'expected',expected
             call die('independent gradient mismatch')
          end if
          checked=checked+1_int64
       end if
       native_velocity(d)%tendency(3*node-2:3*node,k)=expected
    end do
    do d=1,size(native_velocity)
       do s=1,size(native_velocity(d)%gradient)
          if(native_velocity(d)%gradient(s)%active) cycle
          node=native_velocity(d)%gradient(s)%node(1)
          native_velocity(d)%tendency(3*node-2:3*node,k)=0.0_dp
       end do
       native_velocity(d)%ready(k)=.true.
    end do
    pending_layer=0
    write(*,'(a,i0,a,2(1x,i0))') 'Block gradients rank=',rank,' executed/checked=',executed,checked
  end subroutine
end module
