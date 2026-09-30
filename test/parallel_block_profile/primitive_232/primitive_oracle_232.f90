module primitive_oracle_232_mod
  ! Test-only adapter. These Domain reads belong to the independent reference
  ! comparison; this module must never be linked into a production campaign.
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use domain_mod, only: Domain,sol_mean,Float_Field,grid
  use shared_mod, only: S_MASS,N_VARIABLE,zlevels,grav_accel,p_top
  use patch_mod, only: PATCH_SIZE
  use parallel_block_primitive_kernels_mod
  implicit none
  private
  public :: check_regular232,check_southwest232,begin_pentagon232,end_pentagon232
  public :: report232,begin_pressure232,end_pressure232
  public :: check_surface232,check_scalar_diffusion232,check_velocity_diffusion232
  integer(int64) :: counts(10)=0_int64
  real(dp) :: expected_pressure(5)
  real(dp) :: saved_qe(3,0:7),saved_vort(2,0:7)
  integer :: saved_corner=0,saved_ids(0:7)
contains
  subroutine exact(a,b,label)
    real(dp), intent(in) :: a(:),b(:)
    character(*), intent(in) :: label
    if(size(a)/=size(b)) error stop 'primitive oracle extent'
    if(any(transfer(a,[0_int64],size(a))/=transfer(b,[0_int64],size(b)))) then
       write(*,*) 'Primitive232 mismatch: ',label,maxval(abs(a-b))
       error stop 'primitive232 numerical mismatch'
    end if
  end subroutine

  subroutine operands(dom,ids,mass,velocity,k,rho,u,length,pedlen,area,triarea,coriolis)
    type(Domain), intent(in) :: dom
    integer, intent(in) :: ids(0:7),k
    real(dp), intent(in) :: mass(:),velocity(:)
    real(dp), intent(out) :: rho(0:7),u(3,0:7),length(3,0:7),pedlen(3,0:7),area(6,0:7)
    real(dp), intent(out) :: triarea(2,0:7),coriolis(2,0:7)
    integer :: j,id
    do j=0,7
       id=ids(j)
       rho(j)=mass(id+1)+sol_mean(S_MASS,k)%data(dom%id+1)%elts(id+1)
       u(:,j)=velocity(3*id+1:3*id+3)
       length(:,j)=dom%len%elts(3*id+1:3*id+3)
       pedlen(:,j)=dom%pedlen%elts(3*id+1:3*id+3)
       area(:,j)=dom%areas%elts(id+1)%part
       triarea(:,j)=dom%triarea%elts(2*id+1:2*id+2)
       coriolis(:,j)=dom%coriolis%elts(2*id+1:2*id+2)
    end do
  end subroutine

  subroutine check_regular232(dom,ids,mass,velocity,k,context,dual)
    type(Domain), intent(in) :: dom
    integer, intent(in) :: ids(0:7),k
    real(dp), intent(in) :: mass(:),velocity(:),context(2),dual(3)
    real(dp) :: rho(0:7),u(3,0:7),length(3,0:7),pedlen(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
    real(dp) :: qe(3,0:7),vort(2,0:7),ke,b,ctx(2),du(3)
    integer :: id
    call operands(dom,ids,mass,velocity,k,rho,u,length,pedlen,area,triarea,coriolis)
    qe=0.0_dp;vort=0.0_dp;id=ids(0)
    call primitive_regular(rho,u,length,pedlen,area,dom%areas%elts(id+1)%hex_inv,triarea,coriolis, &
         dom%geopot%elts(id+1),dom%geopot_lower%elts(id+1),qe,vort,ke,b,ctx,du)
    call exact(qe(:,0),dom%qe%elts(3*id+1:3*id+3),'regular PV')
    call exact(vort(:,0),dom%vort%elts(2*id+1:2*id+2),'regular vorticity')
    call exact([vort(1,4),vort(2,3)],[dom%vort%elts(2*ids(4)+1),dom%vort%elts(2*ids(3)+2)],'W/S vorticity')
    call exact([ke,b],[dom%ke%elts(id+1),dom%bernoulli%elts(id+1)],'KE/Bernoulli')
    call exact(ctx,context,'PV context');call exact(du,dual,'dual velocity')
    counts(1)=counts(1)+1
  end subroutine

  subroutine check_southwest232(dom,ids,mass,velocity,k,context)
    type(Domain), intent(in) :: dom
    integer, intent(in) :: ids(0:7),k
    real(dp), intent(in) :: mass(:),velocity(:),context(2)
    real(dp) :: rho(0:7),u(3,0:7),length(3,0:7),pedlen(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
    real(dp) :: qe(3,0:7),vort(2,0:7)
    call operands(dom,ids,mass,velocity,k,rho,u,length,pedlen,area,triarea,coriolis)
    qe=0.0_dp;vort=0.0_dp
    call primitive_southwest(rho,u,length,area,triarea,coriolis,context,qe,vort)
    call exact([qe(1,4),qe(2,7),qe(3,3)], &
         [dom%qe%elts(3*ids(4)+1),dom%qe%elts(3*ids(7)+2),dom%qe%elts(3*ids(3)+3)],'SW PV')
    call exact(vort(:,7),dom%vort%elts(2*ids(7)+1:2*ids(7)+2),'SW vorticity')
    counts(2)=counts(2)+1
  end subroutine

  subroutine begin_pentagon232(c,dom,ids,mass,velocity,mean_mass)
    integer, intent(in) :: c,ids(0:7)
    type(Domain), intent(in) :: dom
    real(dp), intent(in) :: mass(:),velocity(:),mean_mass(:)
    real(dp) :: rho(0:7),u(3,0:7),length(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
    integer :: j,id
    if(saved_corner/=0) error stop 'nested primitive pentagon'
    saved_corner=c;saved_ids=ids
    do j=0,7
       id=ids(j)
       rho(j)=mass(id+1)+mean_mass(id+1)
       u(:,j)=velocity(3*id+1:3*id+3)
       length(:,j)=dom%len%elts(3*id+1:3*id+3)
       area(:,j)=dom%areas%elts(id+1)%part
       triarea(:,j)=dom%triarea%elts(2*id+1:2*id+2)
       coriolis(:,j)=dom%coriolis%elts(2*id+1:2*id+2)
       saved_qe(:,j)=dom%qe%elts(3*id+1:3*id+3)
       saved_vort(:,j)=dom%vort%elts(2*id+1:2*id+2)
    end do
    call primitive_pentagon(c,rho,u,length,area,triarea,coriolis,saved_qe,saved_vort)
  end subroutine

  subroutine end_pentagon232(dom)
    type(Domain), intent(in) :: dom
    integer :: qnode(2),qedge(2),vnode(2),vedge(2),j
    select case(saved_corner)
    case(7)
       qnode=[4,3];qedge=[1,3];vnode=[7,7];vedge=[1,2]
    case(6)
       qnode=[0,7];qedge=[1,2];vnode=[7,3];vedge=[1,2]
    case(8)
       qnode=[0,7];qedge=[3,2];vnode=[7,4];vedge=[2,1]
    case(5)
       qnode=[0,0];qedge=[1,3];vnode=[0,0];vedge=[1,2]
    case default
       error stop 'invalid primitive pentagon state'
    end select
    do j=1,2
       call exact([saved_qe(qedge(j),qnode(j))], &
            [dom%qe%elts(3*saved_ids(qnode(j))+qedge(j))],'pentagon PV')
       call exact([saved_vort(vedge(j),vnode(j))], &
            [dom%vort%elts(2*saved_ids(vnode(j))+vedge(j))],'pentagon vorticity')
    end do
    counts(saved_corner-2)=counts(saved_corner-2)+1
    saved_corner=0
  end subroutine

  subroutine begin_pressure232(mass,mean_mass,temp,mean_temp,g,cp,p0,kappa,plower,geopot)
    real(dp), intent(in) :: mass,mean_mass,temp,mean_temp,g,cp,p0,kappa,plower,geopot
    expected_pressure(1)=plower;expected_pressure(2)=geopot
    call pressure_layer(mass,mean_mass,temp,mean_temp,g,cp,p0,kappa, &
         expected_pressure(1),expected_pressure(2),expected_pressure(3),expected_pressure(4),expected_pressure(5))
  end subroutine
  subroutine end_pressure232(plower,geopot,pressure,exner,geopot_lower)
    real(dp), intent(in) :: plower,geopot,pressure,exner,geopot_lower
    call exact(expected_pressure,[plower,geopot,pressure,exner,geopot_lower],'column pressure')
    counts(7)=counts(7)+1
  end subroutine
  subroutine check_surface232(q)
    type(Float_Field), intent(in) :: q(N_VARIABLE,zlevels)
    real(dp) :: m(zlevels),ref(zlevels),pressure
    integer :: d,p,i,k,id
    do d=1,size(grid)
       do p=3,grid(d)%patch%length
          if(grid(d)%patch%elts(p)%deleted) cycle
          do i=1,PATCH_SIZE**2
             id=grid(d)%patch%elts(p)%elts_start+i
             do k=1,zlevels
                m(k)=q(S_MASS,k)%data(d)%elts(id)
                ref(k)=sol_mean(S_MASS,k)%data(d)%elts(id)
             end do
             pressure=surface_pressure(m,ref,grav_accel,p_top)
             call exact([pressure],[grid(d)%surf_press%elts(id)],'surface pressure')
             counts(8)=counts(8)+1
          end do
       end do
    end do
  end subroutine
  subroutine check_scalar_diffusion232(values,length,pedlen,order,nu,expected)
    real(dp), intent(in) :: values(4),length(3),pedlen(3),nu,expected(3)
    integer, intent(in) :: order
    call exact(scalar_diffusion_flux(values,length,pedlen,order,nu),expected,'scalar diffusion')
    counts(9)=counts(9)+1
  end subroutine
  subroutine check_velocity_diffusion232(divergence,vorticity,length,pedlen,div_order,rot_order, &
       div_nu,rot_nu,expected)
    real(dp), intent(in) :: divergence(4),vorticity(4),length(3),pedlen(3),div_nu,rot_nu,expected(3)
    integer, intent(in) :: div_order,rot_order
    call exact(velocity_diffusion_source(divergence,vorticity,length,pedlen,div_order,rot_order, &
         div_nu,rot_nu),expected,'velocity diffusion')
    counts(10)=counts(10)+1
  end subroutine
  subroutine report232
    write(*,'(a,10(1x,i0))') 'Primitive232 exact counts:',counts
  end subroutine
end module
