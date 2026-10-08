  subroutine audit274_begin
    character(len=80) :: filename
    if (audit274_unit/=-1) return
    write(filename,'(a,i0,a)') 'ownership274-rank-',rank,'.txt'
    open(newunit=audit274_unit,file=trim(filename),status='new',action='write')
    write(audit274_unit,*) 'RANK',rank,n_process,N_GLO_DOMAIN,zlevels
    allocate(audit274_samples(N_GLO_DOMAIN))
    audit274_samples=0
  end subroutine

  subroutine audit274_snapshot(d)
    ! Diagnostic-only topology samples: no numerical field reads or mutation.
    integer, intent(in) :: d
    integer :: dom,b,n,i,l,a,s,e,j
    integer, allocatable :: destination(:),block_id(:),node_level(:)
    dom=glo_id(rank+1,d)
    if (dom<0.or.dom>=N_GLO_DOMAIN) error stop 'Stage274 invalid Domain identity'
    if (audit274_samples(dom+1)>=3) return
    audit274_samples(dom+1)=audit274_samples(dom+1)+1
    n=grid(d)%node%length
    allocate(destination(n),block_id(n),node_level(n))
    destination=-1;block_id=-1;node_level=-1
    do b=1,size(audit274_catalog)
       if (audit274_catalog(b)%root_domain/=dom) cycle
       call map_patch(audit274_catalog(b)%root_patch,b)
    end do
    write(audit274_unit,*) 'SNAP',dom,audit274_samples(dom+1),native_velocity(d)%generation, &
         n,level_start,level_end
    do i=1,n
       write(audit274_unit,*) 'NODE',i,destination(i),block_id(i),node_level(i)
    end do
    do l=level_end,level_start,-1
       associate(program=>native_velocity(d)%level(l))
       write(audit274_unit,*) 'LEVEL',l,size(program%action)
       do a=1,size(program%action)
          s=program%operand(a)
          select case(program%action(a))
          case(VELOCITY_DIRECT)
             write(audit274_unit,*) 'DIRECT',a,merge(1,0,program%direct(s)%active),program%direct(s)%node
          case(VELOCITY_RESTRICT)
             write(audit274_unit,*) 'RESTRICT',a,program%restriction(s)%target,program%restriction(s)%child, &
                  program%restriction(s)%neighbor,(merge(1,0,program%restriction(s)%edge(e)),e=1,3)
          case default
             error stop 'Stage274 invalid source action'
          end select
       end do
       end associate
    end do
    write(audit274_unit,*) 'GRADIENTS',size(native_velocity(d)%gradient)
    do j=1,size(native_velocity(d)%gradient)
       write(audit274_unit,*) 'GRADIENT',j,merge(1,0,native_velocity(d)%gradient(j)%active), &
            native_velocity(d)%gradient(j)%node
    end do
    write(audit274_unit,*) 'END'
    flush(audit274_unit)
  contains
    recursive subroutine map_patch(p,catalogue_index)
      integer, intent(in) :: p,catalogue_index
      integer :: first,last,c,child
      if (p<0.or.p>=grid(d)%patch%length) error stop 'Stage274 invalid patch'
      if (grid(d)%patch%elts(p+1)%deleted) return
      first=grid(d)%patch%elts(p+1)%elts_start+1
      last=first+PATCH_SIZE**2-1
      if (first<1.or.last>n) error stop 'Stage274 invalid patch nodes'
      if (any(block_id(first:last)/=-1)) error stop 'Stage274 overlapping patch interiors'
      destination(first:last)=audit274_catalog(catalogue_index)%owner
      block_id(first:last)=audit274_catalog(catalogue_index)%id
      node_level(first:last)=grid(d)%patch%elts(p+1)%level
      do c=1,N_CHDRN
         child=grid(d)%patch%elts(p+1)%children(c)
         if (child>0) call map_patch(child,catalogue_index)
      end do
    end subroutine
  end subroutine
