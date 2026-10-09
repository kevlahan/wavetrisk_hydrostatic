  subroutine prepare_primitive277(generation)
    integer(int64), intent(in) :: generation
    integer :: d,p,i,j,t,a,v,nn,ids(7),offs(N_BDRY+1),dims(2,N_BDRY+1)
    character(len=32) :: setting
    character(len=80) :: filename
    if (.not.p277_initialized) then
       setting='0'
       call get_environment_variable('WAVETRISK_PRIMITIVE277_CHECK',setting)
       if (trim(setting)/='0'.and.trim(setting)/='1') error stop 'Stage277 invalid check selector'
       p277_check=trim(setting)=='1'
       p277_initialized=.true.
       allocate(p277_layer_nodes(zlevels))
       if (rank==0) write(6,'(a,l1)') 'Stage277 cached fused primitive = on; assignment checks = ',p277_check
       if (p277_check) then
          write(filename,'(a,i0,a)') 'primitive277-rank-',rank,'.txt'
          open(newunit=p277_unit,file=trim(filename),status='new',action='write')
          write(p277_unit,*) 'RANK',rank,n_process,N_GLO_DOMAIN,size(grid),zlevels,PATCH_SIZE
       end if
    end if
    p277_transaction=p277_transaction+1
    p277_checks=0;p277_layer_nodes=0
    if (p277_generation==generation) return
    if (allocated(primitive277)) deallocate(primitive277)
    allocate(primitive277(size(grid)))
    do d=1,size(grid)
       allocate(primitive277(d)%patch(grid(d)%patch%length))
       do p=3,grid(d)%patch%length
          if (grid(d)%patch%elts(p)%deleted) cycle
          associate(plan=>primitive277(d)%patch(p))
          allocate(plan%geometry(37,PATCH_SIZE**2),plan%velocity(10,PATCH_SIZE**2),plan%mass(6,PATCH_SIZE**2))
          call get_offs_Domain(grid(d),p-1,offs,dims)
          do j=0,PATCH_SIZE-1
             do i=0,PATCH_SIZE-1
                nn=j*PATCH_SIZE+i+1
                do t=1,7
                   ids(t)=idx(i+P277_DX(t),j+P277_DY(t),offs,dims)
                end do
                plan%mass(:,nn)=ids(1:6)+1
                do t=1,10
                   a=EDGE*ids(P277_VNODE(t))+P277_VEDGE(t)
                   plan%velocity(t,nn)=a
                   plan%geometry(t,nn)=grid(d)%len%elts(a)
                end do
                do t=1,6
                   v=P277_DUAL(t);a=EDGE*ids(P277_VNODE(v))+P277_VEDGE(v)
                   plan%geometry(10+t,nn)=grid(d)%pedlen%elts(a)
                end do
                plan%geometry(17:28,nn)=[grid(d)%areas%elts(ids(1)+1)%part(1),grid(d)%areas%elts(ids(3)+1)%part(3), &
                     grid(d)%areas%elts(ids(6)+1)%part(5),grid(d)%areas%elts(ids(1)+1)%part(2), &
                     grid(d)%areas%elts(ids(6)+1)%part(4),grid(d)%areas%elts(ids(2)+1)%part(6), &
                     grid(d)%areas%elts(ids(5)+1)%part(1),grid(d)%areas%elts(ids(1)+1)%part(3), &
                     grid(d)%areas%elts(ids(2)+1)%part(5),grid(d)%areas%elts(ids(4)+1)%part(2), &
                     grid(d)%areas%elts(ids(3)+1)%part(4),grid(d)%areas%elts(ids(1)+1)%part(6)]
                plan%geometry(29:32,nn)=[grid(d)%coriolis%elts(TRIAG*ids(1)+LORT+1),grid(d)%coriolis%elts(TRIAG*ids(1)+UPLT+1), &
                     grid(d)%coriolis%elts(TRIAG*ids(5)+LORT+1),grid(d)%coriolis%elts(TRIAG*ids(4)+UPLT+1)]
                plan%geometry(33,nn)=grid(d)%areas%elts(ids(1)+1)%hex_inv
                plan%geometry(34:37,nn)=[grid(d)%triarea%elts(TRIAG*ids(1)+LORT+1),grid(d)%triarea%elts(TRIAG*ids(1)+UPLT+1), &
                     grid(d)%triarea%elts(TRIAG*ids(5)+LORT+1),grid(d)%triarea%elts(TRIAG*ids(4)+UPLT+1)]
             end do
          end do
          end associate
       end do
    end do
    p277_generation=generation;p277_rebuilds=p277_rebuilds+1
  end subroutine

  subroutine finish_primitive277
    if (.not.p277_check) return
    if (p277_checks/=19*sum(p277_layer_nodes)) error stop 'Stage277 assignment coverage differs'
    write(p277_unit,*) 'RHS',p277_transaction,p277_generation,p277_rebuilds,p277_checks,p277_layer_nodes
    flush(p277_unit)
  end subroutine

  pure subroutine fused_primitive277(x,g,y)
    ! Direct arithmetic on compact arrays; preserves step1 operation order.
    ! Flux physics is added by the retained adapter after this call.
    real(dp), intent(in) :: x(16),g(37)
    real(dp), intent(out) :: y(19)
    real(dp) :: u(10),v(6),c(4),m(4),pv(4)
    u=x(7:16)*g(1:10)
    v=x(6+P277_DUAL)*g(11:16)
    c(1)=u(1)+u(9)+u(4);c(2)=-(u(4)+u(8)+u(2))
    c(3)=u(3)+u(8)+u(7);c(4)=-(u(1)+u(5)+u(10))
    m(1)=x(1)*g(17)+x(3)*g(18)+x(6)*g(19)
    m(2)=x(1)*g(20)+x(6)*g(21)+x(2)*g(22)
    m(3)=x(5)*g(23)+x(1)*g(24)+x(2)*g(25)
    m(4)=x(4)*g(26)+x(3)*g(27)+x(1)*g(28)
    pv=(g(29:32)+c)/m
    y(1)=v(1)*(0.5_dp*(x(1)+x(3)))
    y(2)=v(2)*(0.5_dp*(x(1)+x(6)))
    y(3)=v(3)*(0.5_dp*(x(1)+x(2)))
    y(4)=0.5_dp*(pv(4)+pv(1));y(5)=0.5_dp*(pv(2)+pv(1));y(6)=0.5_dp*(pv(2)+pv(3))
    y(7)=(u(8)*v(3)+u(4)*v(2)+u(1)*v(1)+u(10)*v(6)+u(6)*v(5)+u(3)*v(4))*g(33)/4
    y(8:11)=c/g(34:37)
    y(12:13)=pv(3:4)
    y(14:19)=v
  end subroutine

  subroutine assert_primitive277(value,reference)
    real(dp), intent(in) :: value(:),reference(:)
    if (size(value)/=size(reference)) error stop 'Stage277 check extent differs'
    if (.not.all(ieee_is_finite(value)).or..not.all(ieee_is_finite(reference))) &
         error stop 'Stage277 nonfinite primitive'
    if (any(transfer(value,[0_int64],size(value))/=transfer(reference,[0_int64],size(reference)))) then
       write(6,*) 'Stage277 rank/RHS',rank,p277_transaction
       write(6,*) 'candidate',value
       write(6,*) 'reference',reference
       error stop 'Stage277 primitive assignment differs bitwise'
    end if
    p277_checks=p277_checks+size(value)
  end subroutine
