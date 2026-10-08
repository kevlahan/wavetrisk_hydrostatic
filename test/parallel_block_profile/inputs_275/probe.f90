  subroutine audit275_transaction_begin
    character(len=80) :: filename
    if (audit275_probe_unit==-1) then
       write(filename,'(a,i0,a)') 'inputs275-rank-',rank,'.txt'
       open(newunit=audit275_probe_unit,file=trim(filename),status='new',action='write')
       write(audit275_probe_unit,*) 'RANK',rank,n_process,N_GLO_DOMAIN,zlevels
    end if
    audit275_transaction=audit275_transaction+1
  end subroutine

  logical function audit275_selected(k) result(selected)
    integer, intent(in) :: k
    selected=(k==1.or.k==zlevels).and.any(audit275_transaction==[1,2,3,9,17,33,41,48])
  end function

  subroutine audit275_source_before(d,l,flux,pv)
    integer, intent(in) :: d,l
    real(dp), intent(in) :: flux(:,:),pv(:,:)
    real(dp), allocatable :: f(:,:),v(:,:),p(:,:)
    integer :: s,j,e,node,n
    n=size(flux,2)
    allocate(f(3,n),v(3,n),p(3,n))
    f=ieee_value(0.0_dp,ieee_quiet_nan);v=f;p=f
    associate(work=>native_velocity(d))
    do s=1,size(work%level(l)%direct)
       if (.not.work%level(l)%direct(s)%active) cycle
       do j=1,19
          e=audit275_component(j);node=work%level(l)%direct(s)%node(audit275_position(j))
          f(e,node)=flux(e,node);v(e,node)=pv(e,node)
       end do
       node=work%level(l)%direct(s)%node(1)
       p(:,node)=work%physics(:,node)
    end do
    audit275_source_shadow=work%source
    call execute_velocity_sources(work%level(l),f,v,p,audit275_source_shadow)
    end associate
  end subroutine

  subroutine audit275_source_after(d)
    integer, intent(in) :: d
    if (any(transfer(audit275_source_shadow,[0_int64],size(audit275_source_shadow))/= &
         transfer(native_velocity(d)%source,[0_int64],size(audit275_source_shadow)))) &
         error stop 'Stage275 minimal source inputs differ'
    audit275_source_checks=audit275_source_checks+1_int64
  end subroutine

  subroutine audit275_gradient_check(d,k,q)
    integer, intent(in) :: d,k
    type(Float_Field), intent(in) :: q(1:N_VARIABLE,1:zlevels)
    real(dp), allocatable :: b(:),x(:),m(:),t(:),result(:,:)
    logical, allocatable :: used(:)
    integer :: s,n,node
    associate(work=>native_velocity(d))
    n=size(work%rho)
    allocate(b(n),x(n),m(n),t(n),used(n),result(3,n))
    b=ieee_value(0.0_dp,ieee_quiet_nan);x=b;m=b;t=b;used=.false.
    do s=1,size(work%gradient)
       if (.not.work%gradient(s)%active) cycle
       used(work%gradient(s)%node)=.true.
    end do
    do node=1,n
       if (.not.used(node)) cycle
       b(node)=grid(d)%bernoulli%elts(node);x(node)=exner_fun(k)%data(d)%elts(node)
       m(node)=work%rho(node);t(node)=work%rho_theta(node)
    end do
    result=reshape(work%tendency(:,k),[3,n])
    call execute_velocity_gradients(work%gradient,b,x,m,t,work%source,result)
    if (any(transfer(result,[0_int64],size(result))/= &
         transfer(work%tendency(:,k),[0_int64],size(result)))) &
         error stop 'Stage275 minimal gradient inputs differ'
    audit275_gradient_checks=audit275_gradient_checks+1_int64
    if (audit275_selected(k)) call audit275_domain_inputs(d,k,q,used)
    end associate
  end subroutine

  subroutine audit275_gradient_begin(k)
    integer, intent(in) :: k
    integer :: context
    if (.not.audit275_selected(k)) return
    context=k
    call audit275_block_consumer(audit275_block_inputs,context)
  end subroutine

  subroutine audit275_gradient_end(k)
    integer, intent(in) :: k
    write(audit275_probe_unit,*) 'CHECK',audit275_transaction,k,audit275_source_checks,audit275_gradient_checks
    flush(audit275_probe_unit)
  end subroutine

  subroutine audit275_block_inputs(catalog_index,block,context)
    integer, intent(in) :: catalog_index
    type(Audit275_Block), intent(in) :: block
    class(*), intent(inout) :: context
    integer :: k,p,i,n,mi,ti,slot,ms,ts
    real(dp) :: values(6)
    select type(context)
    type is(integer)
       k=context
    class default
       error stop 'Stage275 invalid probe context'
    end select
    n=size(block%node);slot=k-block%field_level
    ms=S_MASS-block%scalar_variable;ts=S_TEMP-block%scalar_variable
    if (slot<0.or.slot>=block%n_field_level.or.min(ms,ts)<0.or. &
         max(ms,ts)>=block%n_scalar_variable.or.block%scalar_mult/=1) &
         error stop 'Stage275 unsupported block scalar layout'
    do p=1,size(block%patch)
       if (block%patch(p)%deleted) cycle
       do i=block%patch(p)%elts_start+1,block%patch(p)%elts_start+PATCH_SIZE**2
          mi=(ms*block%n_field_level+slot)*n+i;ti=(ts*block%n_field_level+slot)*n+i
          values=[block%scalar(mi),block%scalar_mean(mi),block%scalar_mean(mi)+block%scalar(mi), &
                  block%scalar(ti),block%scalar_mean(ti),block%scalar_mean(ti)+block%scalar(ti)]
          if (.not.all(ieee_is_finite(values))) error stop 'Stage275 nonfinite block input'
          write(audit275_probe_unit,*) 'B',audit275_transaction,k,block%id,block%patch(p)%level, &
               transfer([block%node(i)%x,block%node(i)%y,block%node(i)%z],[0_int64],3), &
               transfer(values,[0_int64],6)
       end do
    end do
  end subroutine

  subroutine audit275_domain_inputs(d,k,q,used)
    integer, intent(in) :: d,k
    type(Float_Field), intent(in) :: q(1:N_VARIABLE,1:zlevels)
    logical, intent(in) :: used(:)
    integer :: b,dom
    dom=glo_id(rank+1,d)
    do b=1,size(audit275_catalog)
       if (audit275_catalog(b)%root_domain==dom) call visit(audit275_catalog(b)%root_patch,b)
    end do
  contains
    recursive subroutine visit(p,b)
      integer, intent(in) :: p,b
      integer :: i,c,child
      real(dp) :: values(6)
      if (grid(d)%patch%elts(p+1)%deleted) return
      do i=grid(d)%patch%elts(p+1)%elts_start+1,grid(d)%patch%elts(p+1)%elts_start+PATCH_SIZE**2
         if (.not.used(i)) cycle
         values=[q(S_MASS,k)%data(d)%elts(i),sol_mean(S_MASS,k)%data(d)%elts(i),native_velocity(d)%rho(i), &
                 q(S_TEMP,k)%data(d)%elts(i),sol_mean(S_TEMP,k)%data(d)%elts(i),native_velocity(d)%rho_theta(i)]
         if (.not.all(ieee_is_finite(values))) error stop 'Stage275 nonfinite Domain input'
         write(audit275_probe_unit,*) 'D',audit275_transaction,k,audit275_catalog(b)%id, &
              grid(d)%patch%elts(p+1)%level, &
              transfer([grid(d)%node%elts(i)%x,grid(d)%node%elts(i)%y,grid(d)%node%elts(i)%z],[0_int64],3), &
              transfer(values,[0_int64],6)
      end do
      do c=1,N_CHDRN
         child=grid(d)%patch%elts(p+1)%children(c)
         if (child>0) call visit(child,b)
      end do
    end subroutine
  end subroutine
