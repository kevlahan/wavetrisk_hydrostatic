  subroutine p276_begin(k,l)
    integer, intent(in) :: k,l
    integer :: b,d,next,n
    integer(int64) :: generation
    real(dp) :: start
    character(len=80) :: filename
    p276_active=native_velocity_transaction.and.(k==1.or.k==zlevels).and.any(p276_transaction==[1,2,3,9,17,33,41,48])
    if (.not.p276_active) return
    if (p276_unit==-1) then
       write(filename,'(a,i0,a)') 'primitive276-rank-',rank,'.txt'
       open(newunit=p276_unit,file=trim(filename),status='new',action='write')
       write(p276_unit,*) 'RANK',rank,n_process,N_GLO_DOMAIN,size(grid)
    end if
    start=MPI_Wtime()
    p276_prepare=0;p276_capture_time=0;p276_exchange=0;p276_address=0
    p276_gather=0;p276_compute=0;p276_check=0
    p276_used=0;p276_inputs=0;p276_diff=0;p276_outputs=0;p276_output_diff=0
    p276_ref_diff=0;p276_nodes=0;p276_skipped=0
    generation=native_velocity_plan_generation()
    if (generation/=p276_generation) then
       if (allocated(p276_map)) deallocate(p276_map,p276_plan,p276_send,p276_dest)
       allocate(p276_map(size(grid)),p276_plan(size(block_catalog)))
       n=0
       do d=1,size(grid)
          allocate(p276_map(d)%b(grid(d)%patch%length),p276_map(d)%p(grid(d)%patch%length))
          p276_map(d)%b=0;p276_map(d)%p=0
          n=n+grid(d)%patch%length*PATCH_SIZE**2
       end do
       allocate(p276_send(P276_WIDTH,n),p276_dest(n))
       do b=1,size(block_catalog)
          d=loc_id(block_catalog(b)%root_domain+1)+1
          if (d<1.or.d>size(grid)) cycle
          if (glo_id(rank+1,d)/=block_catalog(b)%root_domain) cycle
          next=0
          call map_patch(block_catalog(b)%root_patch)
       end do
       p276_generation=generation
    end if
    p276_prepare=MPI_Wtime()-start
    if (l<level_start.or.l>level_end) error stop 'Stage276 level invalid'
  contains
    recursive subroutine map_patch(p)
      integer, intent(in) :: p
      integer :: c,child
      if (grid(d)%patch%elts(p+1)%deleted) return
      next=next+1
      if (p276_map(d)%b(p+1)/=0) error stop 'Stage276 duplicate patch map'
      p276_map(d)%b(p+1)=b;p276_map(d)%p(p+1)=next
      do c=1,N_CHDRN
         child=grid(d)%patch%elts(p+1)%children(c)
         if (child>0) call map_patch(child)
      end do
    end subroutine
  end subroutine

  subroutine p276_capture(q,d,p,k,physics,reference)
    type(Float_Field), intent(in) :: q(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: d,p,k
    real(dp), intent(in) :: physics(EDGE,PATCH_SIZE**2),reference(9,PATCH_SIZE**2)
    integer :: offs(N_BDRY+1),dims(2,N_BDRY+1),i,j,t,a,ids(7),v,nn,b
    real(dp) :: x(16),g(37),y(9),start
    start=MPI_Wtime()
    b=p276_map(d)%b(p+1)
    if (b==0) then
       p276_skipped=p276_skipped+PATCH_SIZE**2
       return
    end if
    call get_offs_Domain(grid(d),p,offs,dims)
    do j=0,PATCH_SIZE-1
       do i=0,PATCH_SIZE-1
          nn=j*PATCH_SIZE+i+1
          do t=1,7
             ids(t)=idx(i+p276_dx(t),j+p276_dy(t),offs,dims)
          end do
          do t=1,6
             a=ids(t)+1
             x(t)=q(S_MASS,k)%data(d)%elts(a)+sol_mean(S_MASS,k)%data(d)%elts(a)
          end do
          do t=1,10
             a=EDGE*ids(p276_vnode(t))+p276_vedge(t)
             x(6+t)=q(S_VELO,k)%data(d)%elts(a)
             g(t)=grid(d)%len%elts(a)
          end do
          do t=1,6
             v=p276_dual(t);a=EDGE*ids(p276_vnode(v))+p276_vedge(v)
             g(10+t)=grid(d)%pedlen%elts(a)
          end do
          g(17:28)=[grid(d)%areas%elts(ids(1)+1)%part(1),grid(d)%areas%elts(ids(3)+1)%part(3), &
               grid(d)%areas%elts(ids(6)+1)%part(5),grid(d)%areas%elts(ids(1)+1)%part(2), &
               grid(d)%areas%elts(ids(6)+1)%part(4),grid(d)%areas%elts(ids(2)+1)%part(6), &
               grid(d)%areas%elts(ids(5)+1)%part(1),grid(d)%areas%elts(ids(1)+1)%part(3), &
               grid(d)%areas%elts(ids(2)+1)%part(5),grid(d)%areas%elts(ids(4)+1)%part(2), &
               grid(d)%areas%elts(ids(3)+1)%part(4),grid(d)%areas%elts(ids(1)+1)%part(6)]
          g(29:32)=[grid(d)%coriolis%elts(TRIAG*ids(1)+LORT+1),grid(d)%coriolis%elts(TRIAG*ids(1)+UPLT+1), &
               grid(d)%coriolis%elts(TRIAG*ids(5)+LORT+1),grid(d)%coriolis%elts(TRIAG*ids(4)+UPLT+1)]
          g(33)=grid(d)%areas%elts(ids(1)+1)%hex_inv
          g(34:37)=[grid(d)%triarea%elts(TRIAG*ids(1)+LORT+1),grid(d)%triarea%elts(TRIAG*ids(1)+UPLT+1), &
               grid(d)%triarea%elts(TRIAG*ids(5)+LORT+1),grid(d)%triarea%elts(TRIAG*ids(4)+UPLT+1)]
          y=reference(:,nn)
          p276_used=p276_used+1
          p276_dest(p276_used)=block_catalog(b)%owner+1
          p276_send(:,p276_used)=[real([b,p276_map(d)%p(p+1),nn],dp),x,g,physics(:,nn),y]
       end do
    end do
    p276_capture_time=p276_capture_time+MPI_Wtime()-start
  end subroutine

  pure subroutine p276_fused(x,g,physics,y)
    ! Same operation order as step1; no Domain lookup or ownership interpreter.
    real(dp), intent(in) :: x(16),g(37),physics(3)
    real(dp), intent(out) :: y(9)
    real(dp) :: u(10),v(6),c(4),m(4),pv(4)
    u=x(7:16)*g(1:10)
    v=x(6+p276_dual)*g(11:16)
    c(1)=u(1)+u(9)+u(4);c(2)=-(u(4)+u(8)+u(2))
    c(3)=u(3)+u(8)+u(7);c(4)=-(u(1)+u(5)+u(10))
    m(1)=x(1)*g(17)+x(3)*g(18)+x(6)*g(19)
    m(2)=x(1)*g(20)+x(6)*g(21)+x(2)*g(22)
    m(3)=x(5)*g(23)+x(1)*g(24)+x(2)*g(25)
    m(4)=x(4)*g(26)+x(3)*g(27)+x(1)*g(28)
    pv=(g(29:32)+c)/m
    y(1)=v(1)*(0.5_dp*(x(1)+x(3)))+physics(1)
    y(2)=v(2)*(0.5_dp*(x(1)+x(6)))+physics(2)
    y(3)=v(3)*(0.5_dp*(x(1)+x(2)))+physics(3)
    y(4)=0.5_dp*(pv(4)+pv(1));y(5)=0.5_dp*(pv(2)+pv(1));y(6)=0.5_dp*(pv(2)+pv(3))
    y(7)=(u(8)*v(3)+u(4)*v(2)+u(1)*v(1)+u(10)*v(6)+u(6)*v(5)+u(3)*v(4))*g(33)/4
    y(8:9)=c(1:2)/g(34:35)
  end subroutine

  subroutine p276_execute(k,l)
    integer, intent(in) :: k,l
    integer :: sc(n_process),rc(n_process),sd(n_process),rd(n_process),cursor(n_process),i,r,ierr,context
    real(dp), allocatable :: send(:,:)
    real(dp) :: start
    start=MPI_Wtime()
    sc=0
    do i=1,p276_used
       r=p276_dest(i);sc(r)=sc(r)+P276_WIDTH
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) error stop 'Stage276 count exchange failed'
    sd(1)=0;rd(1)=0
    do r=2,n_process
       sd(r)=sd(r-1)+sc(r-1);rd(r)=rd(r-1)+rc(r-1)
    end do
    allocate(send(P276_WIDTH,p276_used))
    if (allocated(p276_recv)) deallocate(p276_recv)
    allocate(p276_recv(P276_WIDTH,sum(rc)/P276_WIDTH))
    cursor=sd/P276_WIDTH
    do i=1,p276_used
       r=p276_dest(i);cursor(r)=cursor(r)+1
       send(:,cursor(r))=p276_send(:,i)
    end do
    call MPI_Alltoallv(send,sc,sd,MPI_DP,p276_recv,rc,rd,MPI_DP,comm,ierr)
    if (ierr/=MPI_SUCCESS) error stop 'Stage276 reference exchange failed'
    p276_exchange=MPI_Wtime()-start
    context=k
    call p276_consumer(p276_block,context)
    if (p276_nodes/=size(p276_recv,2)) error stop 'Stage276 missing block output'
    write(p276_unit,*) 'PHASE',p276_transaction,k,l,p276_generation,p276_used,p276_nodes,p276_skipped, &
         p276_inputs,p276_diff,p276_outputs,p276_output_diff,p276_ref_diff, &
         p276_prepare,p276_capture_time,p276_exchange,p276_address,p276_gather,p276_compute,p276_check
    flush(p276_unit)
  end subroutine

  subroutine p276_block(b,block,context)
    integer, intent(in) :: b
    type(Block276), intent(in) :: block
    class(*), intent(inout) :: context
    integer :: k,p,i,j,t,s,a,n,nn,v,h,row,count_node,ii
    integer, allocatable :: rows(:),nodes(:)
    real(dp), allocatable :: values(:,:),outputs(:,:)
    real(dp) :: start,ref(9),x(16)
    select type(context)
    type is(integer)
       k=context
    class default
       error stop 'Stage276 bad context'
    end select
    start=MPI_Wtime()
    associate(plan=>p276_plan(b))
    n=size(block%patch)*PATCH_SIZE**2
    if (.not.allocated(plan%node)) then
       allocate(plan%storage(7,n),plan%node(7,n),plan%output(9,n))
       do p=1,size(block%patch)
          do j=0,PATCH_SIZE-1
             do i=0,PATCH_SIZE-1
                nn=(p-1)*PATCH_SIZE**2+j*PATCH_SIZE+i+1
                do t=1,7
                   call p276_resolve(block,p,i+p276_dx(t),j+p276_dy(t),plan%storage(t,nn),plan%node(t,nn))
                end do
             end do
          end do
       end do
    end if
    count_node=count(nint(p276_recv(1,:))==b)
    allocate(rows(count_node),nodes(count_node),values(16,count_node),outputs(9,count_node))
    ii=0
    do row=1,size(p276_recv,2)
       if (nint(p276_recv(1,row))/=b) cycle
       ii=ii+1;rows(ii)=row
       p=nint(p276_recv(2,row));nn=nint(p276_recv(3,row))
       if (p<1.or.p>size(block%patch).or.nn<1.or.nn>PATCH_SIZE**2) error stop 'Stage276 packet address invalid'
       nodes(ii)=(p-1)*PATCH_SIZE**2+nn
    end do
    p276_address=p276_address+MPI_Wtime()-start
    start=MPI_Wtime()
    do ii=1,count_node
       nn=nodes(ii)
       do t=1,16
          h=p276_inode(t);v=p276_iedge(t)
          s=plan%storage(h,nn);a=plan%node(h,nn)
          values(t,ii)=p276_value(block,s,a,k,v)
       end do
    end do
    p276_gather=p276_gather+MPI_Wtime()-start
    start=MPI_Wtime()
    do ii=1,count_node
       row=rows(ii)
       call p276_fused(values(:,ii),p276_recv(20:56,row),p276_recv(57:59,row),outputs(:,ii))
       plan%output(:,nodes(ii))=outputs(:,ii)
    end do
    p276_compute=p276_compute+MPI_Wtime()-start
    start=MPI_Wtime()
    do ii=1,count_node
       row=rows(ii);nn=nodes(ii);x=p276_recv(4:19,row)
       call p276_fused(x,p276_recv(20:56,row),p276_recv(57:59,row),ref)
       p276_ref_diff=p276_ref_diff+count(transfer(ref,[0_int64],9)/=transfer(p276_recv(60:68,row),[0_int64],9))
       do t=1,16
          h=p276_inode(t)
          s=plan%storage(h,nn)
          p276_inputs(s)=p276_inputs(s)+1
          if (transfer(values(t,ii),0_int64)/=transfer(x(t),0_int64)) p276_diff(s)=p276_diff(s)+1
       end do
       p276_output_diff=p276_output_diff+count(transfer(outputs(:,ii),[0_int64],9)/= &
            transfer(p276_recv(60:68,row),[0_int64],9))
    end do
    p276_nodes=p276_nodes+count_node;p276_outputs=p276_outputs+9_int64*count_node
    p276_check=p276_check+MPI_Wtime()-start
    end associate
  end subroutine

  subroutine p276_resolve(block,p,i,j,storage,node)
    type(Block276), intent(in) :: block
    integer, intent(in) :: p,i,j
    integer, intent(out) :: storage,node
    integer :: side,record,dims(2),x,y,first,neighbor,link
    if (i>=0.and.i<PATCH_SIZE.and.j>=0.and.j<PATCH_SIZE) then
       storage=STORE_PATCH;node=block%patch(p)%elts_start+j*PATCH_SIZE+i+1
       return
    end if
    if (i<0) then
       side=WEST
       if (j<0) side=SOUTHWEST
       if (j>=PATCH_SIZE) side=NORTHWEST
    else if(i>=PATCH_SIZE) then
       side=EAST
       if (j<0) side=SOUTHEAST
       if (j>=PATCH_SIZE) side=NORTHEAST
    else if(j<0) then
       side=SOUTH
    else
       side=NORTH
    end if
    neighbor=block%patch(p)%neigh(side)
    if (neighbor>0) then
       storage=STORE_PATCH;first=block%patch(neighbor+1)%elts_start;dims=PATCH_SIZE
    else if(neighbor<0) then
       link=-neighbor
       record=block%block_bdry(link)%ghost_id
       if (record>0) then
          storage=STORE_GHOST;first=block%ghost_storage(record)%local_start;dims=PATCH_SIZE
       else
          record=block%block_bdry(link)%storage_id
          storage=STORE_BDRY;first=block%bdry_storage(record)%local_start
          dims=block%block_bdry(link)%dims
       end if
    else
       error stop 'Stage276 missing direct neighbour'
    end if
    x=i;y=j
    if (i<0) x=i+dims(1)
    if (i>=PATCH_SIZE) x=i-PATCH_SIZE
    if (j<0) y=j+dims(2)
    if (j>=PATCH_SIZE) y=j-PATCH_SIZE
    if (x<0.or.x>=dims(1).or.y<0.or.y>=dims(2)) error stop 'Stage276 stencil outside storage'
    node=first+y*dims(1)+x+1
  end subroutine

  real(dp) function p276_value(block,s,a,k,v) result(value)
    type(Block276), intent(in) :: block
    integer, intent(in) :: s,a,k,v
    integer :: n,offset,slot
    slot=k-block%field_level
    if (slot<0.or.slot>=block%n_field_level) error stop 'Stage276 layer outside block'
    select case(s)
    case(STORE_PATCH)
       n=size(block%node)
    case(STORE_BDRY)
       n=size(block%bdry_node)
    case(STORE_GHOST)
       n=size(block%ghost_node)
    case default
       error stop 'Stage276 bad storage'
    end select
    if (a<1.or.a>n) error stop 'Stage276 node outside block'
    if (v==0) then
       offset=((S_MASS-block%scalar_variable)*block%n_field_level+slot)*n+a
       select case(s)
       case(STORE_PATCH)
          value=block%scalar(offset)+block%scalar_mean(offset)
       case(STORE_BDRY)
          value=block%bdry_scalar(offset)+block%bdry_scalar_mean(offset)
       case(STORE_GHOST)
          value=block%ghost_scalar(offset)+block%ghost_scalar_mean(offset)
       end select
    else
       offset=slot*EDGE*n+EDGE*(a-1)+v
       select case(s)
       case(STORE_PATCH)
          value=block%vector(offset)
       case(STORE_BDRY)
          value=block%bdry_vector(offset)
       case(STORE_GHOST)
          value=block%ghost_vector(offset)
       end select
    end if
  end function
