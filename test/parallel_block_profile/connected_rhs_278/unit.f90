program rhs_unit
  use kind_mod, only : dp
  use, intrinsic :: iso_fortran_env, only : int64
  use parallel_block_rhs_mod
  use parallel_block_mass_mod, only : Mass_Level
  use parallel_block_velocity_mod, only : Velocity_Level_Program, Velocity_Gradient_Stencil, &
       VELOCITY_DIRECT,velocity_weights,velocity_gradient,velocity_qperp,velocity_source
  use legacy_primitive_mod, only : legacy_primitive,legacy_pentagon
  implicit none
  type(Rhs_Primitive_Stencil) :: plan(1),repeat_plan(2)
  type(Mass_Level) :: level
  type(Velocity_Level_Program) :: source_plan
  type(Velocity_Gradient_Stencil) :: gradient_plan(1)
  real(dp) :: m(9,3),t(9,3),mm(9,3),mt(9,3),u(12),g(37),x(16),reference(11)
  real(dp) :: p(9,3),ex(9,3),phi(9,3),surface(9),physics(3,9),flux(3,9),pv(3,9),ke(9),b(9),vort(6)
  real(dp) :: sw_x(3),sw_g(12),sw_y(8),sw_physics(3,1)
  real(dp) :: flat_flux(27),rhs(9),velocity_rhs(3,9),stage_m(9),stage_u(3,9),base_u(3,9),source(3,9)
  real(dp) :: lower,upper,geopot,old_geopot,rho,theta,refp,refex,refphi,parts(6,4),inv(4)
  real(dp) :: expected_source(3),expected_velocity(3),diff(3),expected_diff(3)
  integer :: i,k,trial,order
  character(len=32) :: mode
  call get_command_argument(1,mode)
  plan(1)%target=1;plan(1)%mass=[1,2,3,4,5,6];plan(1)%velocity=[1,2,3,4,5,6,7,8,9,10]
  plan(1)%vorticity=[1,2,3,4]
  allocate(level%restriction(0),level%divergence(1))
  level%divergence(1)%target=1;level%divergence(1)%active=.true.
  level%divergence(1)%edge=[1,13,2,20,3,12];level%divergence(1)%inverse_area=0.03125_dp
  allocate(source_plan%action(1),source_plan%operand(1),source_plan%direct(1),source_plan%restriction(0))
  source_plan%action=VELOCITY_DIRECT;source_plan%operand=1
  source_plan%direct(1)%node=[1,2,3,4,5,6,7,8,9];source_plan%direct(1)%active=.true.
  gradient_plan(1)%node=[1,3,6,2];gradient_plan(1)%active=.true.
  if (trim(mode)=='bad_column') then
     call rhs_columns(m,t,mm,mt,surface(1:8),9.81_dp,1004.0_dp,100000.0_dp,0.2857_dp,1000.0_dp,p,ex,phi)
     error stop 'bad extent accepted'
  end if
  if (trim(mode)=='bad_rk') then
     call rhs_rk(m(:,1),base_u(:,1:8),rhs,velocity_rhs(:,1:8),0.125_dp,stage_m,stage_u(:,1:8))
     error stop 'bad RK layout accepted'
  end if
  if (trim(mode)=='missing_sw') then
     plan(1)%southwest=.true.
     call rhs_primitives(plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort)
     error stop 'missing SW physics accepted'
  end if
  if (trim(mode)=='bad_sw') then
     plan(1)%southwest=.true.;plan(1)%sw_mass=7;plan(1)%sw_velocity=[11,12]
     plan(1)%sw_edge=reshape([1,5,2,7,3,0],[2,3]);plan(1)%sw_vorticity=[5,6]
     call rhs_primitives(plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort,sw_physics)
     error stop 'bad SW stencil accepted'
  end if
  if (trim(mode)=='bad_pentagon') call test_pentagons(0)
  do trial=1,512
     call test_pentagons(trial)
     call test_pressure_visits(trial)
     do k=1,3
        do i=1,9
           m(i,k)=13.0_dp+sin(real(trial+7*i+11*k,dp))
           mm(i,k)=2000.0_dp+17*i+k
           t(i,k)=300*m(i,k)+cos(real(trial+i+k,dp))
           mt(i,k)=300*mm(i,k)+i
        end do
     end do
     do i=1,12
        u(i)=sin(real(17*trial+i,dp))*3
     end do
     do i=1,37
        g(i)=0.5_dp+abs(cos(real(trial+5*i,dp)))
     end do
     surface=100.0_dp
     call rhs_columns(m,t,mm,mt,surface,9.81_dp,1004.0_dp,100000.0_dp,0.2857_dp,1000.0_dp,p,ex,phi)
     do i=1,9
        lower=0.0_dp
        do k=1,3
           rho=mm(i,k)+m(i,k);lower=lower+rho
        end do
        lower=9.81_dp*lower+1000.0_dp;geopot=surface(i)
        do k=1,3
           rho=mm(i,k)+m(i,k);theta=mt(i,k)+t(i,k)
           upper=lower-9.81_dp*rho;refp=0.5_dp*(lower+upper)
           refex=1004.0_dp*(refp/100000.0_dp)**0.2857_dp
           old_geopot=geopot;geopot=old_geopot+9.81_dp*0.2857_dp*theta*refex/refp
           refphi=0.5_dp*(geopot+old_geopot)
           call exact([p(i,k),ex(i,k),phi(i,k)],[refp,refex,refphi],'column')
           lower=upper
        end do
     end do
     plan(1)%geometry=g
     physics=0.0625_dp;flux=0.0_dp;pv=0.0_dp;ke=0.0_dp;b=0.0_dp;vort=0.0_dp
     if (trim(mode)=='bad_stencil') plan(1)%velocity(1)=0
     call rhs_primitives(plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort)
     x(1:6)=m(1:6,1)+mm(1:6,1);x(7:16)=u(1:10)
     call legacy_primitive(x,g,reference)
     call exact([flux(:,1),pv(:,1),ke(1),vort(1:4)],[reference(1:3)+physics(:,1),reference(4:11)],'primitive')
     call exact([b(1)],[reference(7)+phi(1,1)],'Bernoulli')
     ! comp_SW must reuse the immediately preceding comput's corner PV and
     ! dual velocities. Its three edge destinations can live in distinct slots.
     plan(1)%southwest=.true.;plan(1)%sw_mass=7;plan(1)%sw_velocity=[11,12]
     plan(1)%sw_edge=reshape([1,5,2,7,3,4],[2,3]);plan(1)%sw_vorticity=[5,6]
     sw_g=g([1,2,17,18,19,20,21,22,29,30,34,35])*1.25_dp
     plan(1)%sw_geometry=sw_g
     sw_x=[m(7,1)+mm(7,1),u(11:12)]
     sw_physics(:,1)=[0.000125_dp,-0.00025_dp,0.000375_dp]
     call legacy_primitive(x,g,reference,sw_x,sw_g,sw_y)
     call rhs_primitives(plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort,sw_physics)
     call exact([flux(1,5),flux(2,7),flux(3,4)],sw_y(1:3)+sw_physics(:,1),'SW flux')
     call exact([pv(1,5),pv(2,7),pv(3,4)],sw_y(4:6),'SW PV')
     call exact(vort(5:6),sw_y(7:8),'SW vorticity')
     ! Overlapping output aliases: SW follows direct, with exact final writes.
     plan(1)%sw_edge=reshape([1,1,2,1,3,1],[2,3]);plan(1)%sw_vorticity=[1,2]
     call rhs_primitives(plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort,sw_physics)
     call exact(flux(:,1),sw_y(1:3)+sw_physics(:,1),'SW overwrites direct flux')
     call exact(pv(:,1),sw_y(4:6),'SW overwrites direct PV')
     call exact(vort(1:2),sw_y(7:8),'SW final triangle writer')
     ! Disabled SW must never read its sentinel addresses or physics input.
     plan(1)%southwest=.false.;plan(1)%sw_mass=0;plan(1)%sw_velocity=0;plan(1)%sw_edge=0
     ! Consecutive producers overwrite shared triangle slots; final writer wins.
     repeat_plan(1)=plan(1);repeat_plan(2)=plan(1)
     repeat_plan(2)%velocity=[10,9,8,7,6,5,4,3,2,1]
     call rhs_primitives(repeat_plan,m(:,1),mm(:,1),u,phi(:,1),physics,flux,pv,ke,b,vort)
     x(7:16)=u(10:1:-1);call legacy_primitive(x,g,reference)
     call exact(vort(1:4),reference(8:11),'final shared vorticity')
     ! Restriction is ordered injection. The second read sees the first write.
     call rhs_scalar_restrict([2,3],[1,2],b,ex(:,1))
     call exact([b(2),b(3)],[b(1),b(1)],'Bernoulli injection')
     call exact([ex(2,1),ex(3,1)],[ex(1,1),ex(1,1)],'Exner injection')
     do i=1,9
        flux(:,i)=flux(:,1)+real(i,dp)*0.125_dp
        pv(:,i)=pv(:,1)+real(i,dp)*0.000125_dp
     end do
     rhs=0.0_dp
     flat_flux=reshape(flux,[27])
     call rhs_mass_restrict(level,flat_flux,rhs)
     call rhs_mass_divergence(level,flat_flux,rhs)
     call exact([rhs(1)],[-(flux(1,1)-flux(1,5)+flux(2,1)-flux(2,7)+flux(3,1)-flux(3,4))*0.03125_dp], &
          'mass divergence')
     parts=0.0625_dp;inv=1.0_dp
     source_plan%direct(1)%weights=velocity_weights(parts,inv)
     source_plan%direct(1)%length=[2.0_dp,3.0_dp,4.0_dp]
     source=0.0_dp
     call rhs_velocity_source(source_plan,flux,pv,physics,source)
     expected_source=velocity_source(velocity_qperp(flux,pv,source_plan%direct(1)%weights), &
          physics(:,1),source_plan%direct(1)%length,.true.)
     call exact(source(:,1),expected_source,'connected velocity source')
     gradient_plan(1)%length=source_plan%direct(1)%length
     velocity_rhs=0.0_dp
     call rhs_velocity_gradient(gradient_plan,b,ex(:,1),m(:,1)+mm(:,1),t(:,1)+mt(:,1),source,velocity_rhs)
     expected_velocity=velocity_gradient(expected_source,gradient_plan(1)%length,b([1,3,6,2]),ex([1,3,6,2],1), &
          m([1,3,6,2],1)+mm([1,3,6,2],1),t([1,3,6,2],1)+mt([1,3,6,2],1),.true.)
     call exact(velocity_rhs(:,1),expected_velocity,'connected gradient')
     base_u=1.25_dp
     call rhs_rk(m(:,1),base_u,rhs,velocity_rhs,0.125_dp,stage_m,stage_u)
     call exact(stage_m,m(:,1)+0.125_dp*rhs,'mass RK')
     call exact(stage_u(:,1),base_u(:,1)+0.125_dp*expected_velocity,'velocity RK')
     do order=0,2
        diff=rhs_scalar_diffusion(m(1:4,1),[2.0_dp,3.0_dp,4.0_dp],[5.0_dp,6.0_dp,7.0_dp],0.25_dp,order)
        expected_diff=0.0_dp
        if (order/=0) expected_diff=(-1)**order*0.25_dp* &
             [(m(2,1)-m(1,1))/2,(m(1,1)-m(3,1))/3,(m(4,1)-m(1,1))/4]*[5.0_dp,6.0_dp,7.0_dp]
        call exact(diff,expected_diff,'scalar diffusion')
        diff=rhs_velocity_diffusion(m(1:4,1),vort(1:4),[2.0_dp,3.0_dp,4.0_dp],[5.0_dp,6.0_dp,7.0_dp], &
             0.25_dp,0.125_dp,order,order)
        expected_diff=0.0_dp
        if (order/=0) then
           expected_diff=(-1)**(order-1)*0.25_dp* &
                [(m(2,1)-m(1,1))/2,(m(1,1)-m(3,1))/3,(m(4,1)-m(1,1))/4]
           expected_diff=expected_diff-(-1)**(order-1)*0.125_dp* &
                [(vort(1)-vort(4))/5,(vort(1)-vort(2))/6,(vort(3)-vort(2))/7]
        end if
        call exact(diff,expected_diff,'velocity diffusion')
     end do
  end do
  print *, 'PASS: 512 connected trials; direct/SW and all four pentagon oracles, aliases, pressure, restriction, RHS and RK'
contains
  subroutine test_pressure_visits(trial)
    integer,intent(in)::trial
    integer,parameter::visits(4)=[2,1,0,3],nodes(6)=[1,4,2,1,4,4]
    real(dp)::m(4,3),t(4,3),mm(4,3),mt(4,3),surface(4),lower(4),phi(4),p(4),ex(4),mid(4)
    real(dp)::ref_lower(4),ref_phi(4),ref_p(4),ref_ex(4),ref_mid(4),upper,old,rho,theta
    integer::n,k,j
    do k=1,3
       do n=1,4
          m(n,k)=sin(real(trial+3*n+7*k,dp));mm(n,k)=100.0_dp+13*n+3*k
          t(n,k)=cos(real(trial+5*n+11*k,dp));mt(n,k)=10000.0_dp+19*n+7*k
       end do
    end do
    surface=[3.0_dp,7.0_dp,11.0_dp,13.0_dp];phi=-99.0_dp
    call rhs_column_start(visits,m,mm,surface,9.81_dp,10000.0_dp,lower,phi)
    ref_lower=0.0_dp;ref_phi=-99.0_dp
    do k=1,3
       do j=1,size(nodes)
          n=nodes(j);ref_lower(n)=ref_lower(n)+(mm(n,k)+m(n,k))
          ref_phi(n)=surface(n)
       end do
    end do
    ref_lower=9.81_dp*ref_lower+10000.0_dp
    call exact(lower,ref_lower,'surface repeated visits')
    call exact(phi,ref_phi,'surface untouched alias')
    p=-99.0_dp;ex=-99.0_dp;mid=-99.0_dp;ref_p=p;ref_ex=ex;ref_mid=mid
    do k=1,3
       ! Separate patches must retain their intermediate lower/upper state.
       call rhs_pressure_patch(nodes(:3),m(:,k),t(:,k),mm(:,k),mt(:,k), &
            9.81_dp,1004.0_dp,100000.0_dp,0.2857_dp,lower,phi,p,ex,mid)
       call rhs_pressure_patch(nodes(4:),m(:,k),t(:,k),mm(:,k),mt(:,k), &
            9.81_dp,1004.0_dp,100000.0_dp,0.2857_dp,lower,phi,p,ex,mid)
       do j=1,size(nodes)
          n=nodes(j);rho=mm(n,k)+m(n,k);theta=mt(n,k)+t(n,k)
          upper=ref_lower(n)-9.81_dp*rho
          ref_p(n)=0.5_dp*(ref_lower(n)+upper)
          ref_ex(n)=1004.0_dp*(ref_p(n)/100000.0_dp)**0.2857_dp
          old=ref_phi(n)
          ref_phi(n)=old+9.81_dp*0.2857_dp*theta*ref_ex(n)/ref_p(n)
          ref_mid(n)=0.5_dp*(ref_phi(n)+old);ref_lower(n)=upper
       end do
       call exact(lower,ref_lower,'pressure repeated visits lower')
       call exact(phi,ref_phi,'pressure repeated visits geopotential')
       call exact(p,ref_p,'pressure repeated visits pressure')
       call exact(ex,ref_ex,'pressure repeated visits exner')
       call exact(mid,ref_mid,'pressure repeated visits midpoint')
    end do
  end subroutine
  subroutine test_pentagons(trial)
    integer, intent(in) :: trial
    type(Rhs_Pentagon_Stencil) :: plan(4)
    real(dp) :: mass(7),mean_mass(7),velocity(3,7),vector(21),length(3,7),parts(6,7),area(2,7),coriolis(2,7)
    real(dp) :: pv(3,7),vort(14),reference_pv(21),reference_vort(14),initial_pv(3,7),initial_vort(14)
    integer :: n,e,c
    do n=1,7
       mass(n)=1.25_dp+sin(real(trial+3*n,dp));mean_mass(n)=1000.0_dp+17*n
       do e=1,3
          velocity(e,n)=sin(real(trial+3*n+7*e,dp))
          length(e,n)=0.5_dp+abs(cos(real(trial+5*n+11*e,dp)))
          initial_pv(e,n)=sin(real(trial+13*n+3*e,dp))*0.0001_dp
       end do
       do e=1,6
          parts(e,n)=0.01_dp+abs(sin(real(trial+7*n+19*e,dp)))
       end do
       do e=1,2
          area(e,n)=0.5_dp+abs(sin(real(trial+11*n+17*e,dp)))
          coriolis(e,n)=cos(real(trial+13*n+5*e,dp))*0.0001_dp
          initial_vort(2*(n-1)+e)=sin(real(trial+23*n+11*e,dp))*0.0002_dp
       end do
    end do
    vector=reshape(velocity,[21])
    do c=1,4
       plan(c)%corner=c+4;plan(c)%mass=[1,2,3,4,5,6,7];plan(c)%pv_node=[1,2,3,4,5,6,7]
       plan(c)%velocity=reshape([(n,n=1,21)],[3,7])
       plan(c)%vorticity=reshape([(n,n=1,14)],[2,7])
       plan(c)%length=length;plan(c)%parts=parts;plan(c)%triangle_area=area;plan(c)%coriolis=coriolis
       pv=initial_pv;vort=initial_vort;reference_pv=reshape(pv,[21]);reference_vort=vort
       if (trial==0) plan(c)%corner=4
       call rhs_pentagons(plan(c:c),mass,mean_mass,vector,pv,vort)
       call legacy_pentagon(c+4,mass+mean_mass,velocity,length,parts,area,coriolis,reference_pv,reference_vort)
       call exact(reshape(pv,[21]),reference_pv,'pentagon PV full field')
       call exact(vort,reference_vort,'pentagon vorticity full field')
    end do
    ! Multiple corrections read previous corrections' final shared values.
    pv=initial_pv;vort=initial_vort;reference_pv=reshape(pv,[21]);reference_vort=vort
    call rhs_pentagons(plan,mass,mean_mass,vector,pv,vort)
    do c=1,4
       call legacy_pentagon(c+4,mass+mean_mass,velocity,length,parts,area,coriolis,reference_pv,reference_vort)
    end do
    call exact(reshape(pv,[21]),reference_pv,'ordered pentagon PV')
    call exact(vort,reference_vort,'ordered pentagon vorticity')
  end subroutine test_pentagons

  subroutine exact(a,b,label)
    real(dp), intent(in) :: a(:),b(:)
    character(len=*), intent(in) :: label
    if (size(a)/=size(b)) error stop 'comparison extent'
    if (any(transfer(a,[0_int64],size(a))/=transfer(b,[0_int64],size(b)))) then
       print *, label, a, b
       error stop 'bit-for-bit comparison failed'
    end if
  end subroutine
end program
