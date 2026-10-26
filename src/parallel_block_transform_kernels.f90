module parallel_block_transform_kernels_mod
  ! Fixed compact numerical kernels. Geometry compilation is separate from
  ! execution; runtime loops never inspect a Domain or resolve an owner.
  use iso_fortran_env,only : int64
  use kind_mod,only : dp
  use shared_mod,only : EDGE,N_BDRY,N_CHDRN,level_start,level_end,FROZEN,ADJZONE, &
       zlevels,DG,RT,UP,LORT,UPLT,TRIAG,NORTHEAST,NORTHWEST,IMINUSJPLUS,IPLUSJMINUS,IJMINUS, &
       end_pt,opp_no,hex_s_offs,hex_sides
  use domain_mod,only : grid,idx,ed_idx,get_offs_Domain,chd_offs
  use parallel_block_rhs_placement_mod,only : compile_rhs_owner_map
  use arch_mod,only : rank,glo_id
  use patch_mod,only : PATCH_SIZE,LAST
  implicit none
  private
  public :: check_forward_transform
  public :: Transform_Scalar,Transform_Inner,Transform_Outer,Transform_Geometry
  public :: compile_transform_geometry,execute_transform_scalar,execute_transform_inner,execute_transform_outer
  type :: Transform_Scalar
     integer :: d=0,global_domain=0,level=0,phase=0,target=0,base=0,area=0,mask=0,n=0
     integer :: source(11)=0,geometry(11)=0,inverse_active=0,forward_active=0
     real(dp) :: weight(11)=0,inverse_area=0
  end type
  type :: Transform_Inner
     integer :: d=0,global_domain=0,level=0,cell=0,node(32)=0,edge(32)=0,source(32)=0,target(6)=0
     integer :: forward_active(6)=0
     real(dp) :: geometry(32)=0
  end type
  type :: Transform_Outer
     integer :: global_domain=0,level=0,kind=0,cell=0,component=0,source(13)=0,target(3)=0,forward_active=0
     real(dp) :: weight(9)=0
  end type
  type :: Transform_Geometry
     type(Transform_Scalar),allocatable :: scalar(:)
     type(Transform_Inner),allocatable :: inner(:)
     type(Transform_Outer),allocatable :: outer(:)
  end type
  real(dp),parameter :: base_weight(9)= &
       [16.0_dp,-1.0_dp,1.0_dp,1.0_dp,-1.0_dp,-1.0_dp,-1.0_dp,1.0_dp,1.0_dp]/16.0_dp
contains
  subroutine die(message)
    character(*),intent(in) :: message
    write(*,*) 'transform geometry: ',message
    error stop 'transform geometry invalid'
  end subroutine
  integer function node_slot(d,id) result(slot)
    integer,intent(in) :: d,id
    if(d<1.or.d>size(grid)) call die('Domain index')
    if(id<0.or.id>=grid(d)%node%length) call die('node index')
    slot=id+1
  end function
  integer function edge_slot(d,id) result(slot)
    integer,intent(in) :: d,id
    if(d<1.or.d>size(grid)) call die('Domain index')
    if(id<0.or.id>=EDGE*grid(d)%node%length) call die('edge index')
    slot=id+1
  end function
  subroutine compile_transform_geometry(program)
    type(Transform_Geometry),intent(out) :: program
    integer(int64) :: old_bytes,new_bytes,total
    integer :: d
    character(8) :: option
    call compile_transform_cells(program)
    call compile_transform_outer(program)
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
    if(trim(option)=='1') then
       total=sum([(int(grid(d)%node%length,int64),d=1,size(grid))])
       old_bytes=total*(2_int64*storage_size(Transform_Scalar())+storage_size(Transform_Inner())+ &
            EDGE*int(storage_size(Transform_Outer()),int64))/8
       new_bytes=(size(program%scalar,kind=int64)*storage_size(Transform_Scalar())+ &
            size(program%inner,kind=int64)*storage_size(Transform_Inner())+ &
            size(program%outer,kind=int64)*storage_size(Transform_Outer()))/8
       write(*,*) 'Structure284 tapes:',rank,old_bytes,new_bytes,size(program%scalar),size(program%inner), &
            size(program%outer)
    end if
    call dump_transform_dependencies(program)
  end subroutine
  subroutine dump_transform_dependencies(program)
    type(Transform_Geometry),intent(in) :: program
    integer,allocatable :: owners(:),keys(:,:)
    integer :: d,u,i
    character(100) :: filename
    character(8) :: setting
    call get_environment_variable('WAVETRISK_TRANSFORM_KERNEL_CHECK',setting)
    if(trim(setting)/='1') return
    do d=1,size(grid)
       call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),owners,keys)
       write(filename,'(a,i0,a)') 'transform-dependencies-',glo_id(rank+1,d),'.txt'
       open(newunit=u,file=trim(filename),status='replace')
       write(u,*) size(owners),level_start,level_end
       do i=1,size(owners)
          write(u,*) 'N',i,owners(i),keys(:,i)
       end do
       do i=1,size(program%scalar)
          associate(op=>program%scalar(i))
          if(op%d/=d) cycle
          write(u,*) 'S',op%level,op%phase,op%target,op%base,op%source,op%inverse_active,op%forward_active
          end associate
       end do
       do i=1,size(program%inner)
          associate(op=>program%inner(i))
          if(op%d/=d) cycle
          write(u,*) 'I',op%level,op%target,op%source
          end associate
       end do
       do i=1,size(program%outer)
          associate(op=>program%outer(i))
          if(op%global_domain/=glo_id(rank+1,d)+1) cycle
          write(u,*) 'O',op%level,op%kind,op%target,op%source,op%forward_active
          end associate
       end do
       close(u)
    end do
  end subroutine

  subroutine compile_transform_cells(program)
    type(Transform_Geometry),intent(inout) :: program
    integer :: total,ncell
    integer :: d,l,p,pi,c,ch,i,j,x,y,ip,jp,n,m,t,a,offp(N_BDRY+1),offc(N_BDRY+1)
    integer :: dimp(2,N_BDRY+1),dimc(2,N_BDRY+1),sx(11),sy(11),q,baseid
    integer :: tx(4),ty(4),wi(11)
    type(Transform_Scalar) :: op
    type(Transform_Inner) :: cell
    total=sum([(grid(d)%node%length,d=1,size(grid))]);ncell=0
    do l=level_start,level_end-1
       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             ncell=ncell+count(grid(d)%patch%elts(p+1)%children>0)*(PATCH_SIZE/2)**2
          end do
       end do
    end do
    if (allocated(program%scalar)) deallocate(program%scalar,program%inner)
    allocate(program%scalar(4*ncell),program%inner(ncell))
    n=0;m=0
    wi=[1,2,3,1,2,3,4,1,2,3,4]
    do l=level_start,level_end-1
       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             if (.not.any(grid(d)%patch%elts(p+1)%children>0)) cycle
             call get_offs_Domain(grid(d),p,offp,dimp)
             do c=1,N_CHDRN
                ch=grid(d)%patch%elts(p+1)%children(c)
                if(ch==0) cycle
                call get_offs_Domain(grid(d),ch,offc,dimc)
                do j=1,PATCH_SIZE/2
                   y=2*(j-1);jp=j-1+chd_offs(2,c)
                   do i=1,PATCH_SIZE/2
                      x=2*(i-1);ip=i-1+chd_offs(1,c)
                      tx=[x,x+1,x+1,x];ty=[y,y,y+1,y+1]
                      do t=1,4
                         op=Transform_Scalar(d=d,level=l,phase=merge(1,2,t==1))
                         op%target=node_slot(d,idx(tx(t),ty(t),offc,dimc))
                         op%mask=idx(x,y,offc,dimc)
                         if(t==1) then
                            op%n=11;sx=x+[1,1,2,0,-1,-1,-2,-1,0,-1,1]
                            sy=y+[0,1,1,1,0,1,-1,-1,-1,-2,-1]
                            baseid=idx(ip,jp,offp,dimp)
                            op%base=node_slot(d,baseid);op%area=baseid
                         else
                            op%n=4
                            select case(t)
                            case(2)
                               sx(:4)=[x,x+2,x+2,x];sy(:4)=[y,y,y+2,y-2]
                            case(3)
                               sx(:4)=[x+2,x,x+2,x];sy(:4)=[y+2,y,y,y+2]
                            case(4)
                               sx(:4)=[x,x,x-2,x+2];sy(:4)=[y,y+2,y,y+2]
                            end select
                            op%base=op%target;op%area=idx(tx(t),ty(t),offc,dimc)
                         end if
                         do q=1,op%n
                            a=idx(sx(q),sy(q),offc,dimc)
                            op%source(q)=node_slot(d,a)
                            op%geometry(q)=merge(a,op%area,t==1)
                         end do
                         op%global_domain=glo_id(rank+1,d)+1
                         op%inverse_active=merge(1,0,t/=1.or.grid(d)%mask_n%elts(op%mask+1)/=FROZEN)
                         op%forward_active=merge(1,0,t/=1.and.grid(d)%mask_n%elts(op%area+1)>=ADJZONE)
                         op%inverse_area=grid(d)%areas%elts(op%area+1)%hex_inv
                         do q=1,op%n
                            op%weight(q)=grid(d)%overl_areas%elts(op%geometry(q)+1)%a(merge(wi(q),q,t==1))
                         end do
                         n=n+1
                         if(n>size(program%scalar)) call die('scalar tape capacity')
                         program%scalar(n)=op
                      end do
                      cell=Transform_Inner(d=d,level=l)
                      cell%cell=idx(x,y,offc,dimc)+1
                      do t=LORT,UPLT
                         call edge_geom(4*(t-LORT)+1,ip,jp,DG,.true.)
                         call edge_geom(4*(t-LORT)+2,ip-t+1,jp,UP,.true.)
                         call edge_geom(4*(t-LORT)+3,ip,jp+t,RT,.true.)
                         call triangle_geom(4*(t-LORT)+4,ip,jp,t,.true.)
                      end do
                      call triangle_geom(9,x,y,LORT,.false.)
                      call edge_geom(10,x,y,RT,.false.)
                      call edge_geom(11,x,y,DG,.false.)
                      call edge_geom(12,x+1,y,UP,.false.)
                      call triangle_geom(13,x+1,y,LORT,.false.)
                      call edge_geom(14,x+1,y,RT,.false.)
                      call edge_geom(15,x+2,y,UP,.false.)
                      call edge_geom(16,x+1,y,DG,.false.)
                      call triangle_geom(17,x+1,y+1,LORT,.false.)
                      call edge_geom(18,x+1,y+1,DG,.false.)
                      call edge_geom(19,x+2,y+1,UP,.false.)
                      call edge_geom(20,x+1,y+1,RT,.false.)
                      call triangle_geom(21,x,y,UPLT,.false.)
                      call edge_geom(22,x,y,UP,.false.)
                      call edge_geom(23,x,y,DG,.false.)
                      call edge_geom(24,x,y+1,RT,.false.)
                      call triangle_geom(25,x,y+1,UPLT,.false.)
                      call edge_geom(26,x,y+1,UP,.false.)
                      call edge_geom(27,x,y+2,RT,.false.)
                      call edge_geom(28,x,y+1,DG,.false.)
                      call triangle_geom(29,x+1,y+1,UPLT,.false.)
                      call edge_geom(30,x+1,y+1,DG,.false.)
                      call edge_geom(31,x+1,y+2,RT,.false.)
                      call edge_geom(32,x+1,y+1,UP,.false.)
                      cell%target=cell%source([12,16,20,24,28,32])
                      cell%global_domain=glo_id(rank+1,d)+1
                      do q=1,32
                         if(any([4,8,9,13,17,21,25,29]==q)) then
                            cell%geometry(q)=grid(d)%triarea%elts(TRIAG*cell%node(q)+cell%edge(q)+1)
                         else
                            cell%geometry(q)=grid(d)%len%elts(EDGE*cell%node(q)+cell%edge(q)+1)
                         end if
                      end do
                      do q=1,6
                         cell%forward_active(q)=merge(1,0,grid(d)%mask_e%elts(cell%target(q))>=ADJZONE)
                      end do
                      m=m+1
                      if(m>size(program%inner)) call die('inner tape capacity')
                      program%inner(m)=cell
                   end do
                end do
             end do
          end do
       end do
    end do
    if(n/=4*ncell.or.m/=ncell) call die('cell tape count')
  contains
    subroutine edge_geom(slot,xx,yy,e,parent)
      integer, intent(in) :: slot,xx,yy,e
      logical, intent(in) :: parent
      integer :: id
      if(parent) then
         id=idx(xx,yy,offp,dimp)
      else
         id=idx(xx,yy,offc,dimc)
      end if
      cell%node(slot)=id;cell%edge(slot)=e
      cell%source(slot)=edge_slot(d,EDGE*id+e)
    end subroutine
    subroutine triangle_geom(slot,xx,yy,e,parent)
      integer, intent(in) :: slot,xx,yy,e
      logical, intent(in) :: parent
      if(parent) then
         cell%node(slot)=idx(xx,yy,offp,dimp)
      else
         cell%node(slot)=idx(xx,yy,offc,dimc)
      end if
      cell%edge(slot)=e
    end subroutine
  end subroutine

  subroutine compile_transform_outer(program)
    type(Transform_Geometry),intent(inout) :: program
    integer :: noperation,total,capacity
    integer :: l,d,p,pi,c,child,i,j,ic,jc,ip,jp,e,corner,q,a,b,t,mode
    integer :: op(N_BDRY+1),oc(N_BDRY+1),dp0(2,N_BDRY+1),dc(2,N_BDRY+1)
    type(Transform_Outer) :: item
    total=sum([(grid(d)%node%length,d=1,size(grid))]);noperation=0
    capacity=0
    do l=level_start-1,level_end-1
       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             capacity=capacity+count(grid(d)%patch%elts(p+1)%children>0)*(UP-RT+1)*(PATCH_SIZE/2+1)**2
          end do
          do corner=NORTHEAST,NORTHWEST
             if(.not.grid(d)%penta(corner)) cycle
             p=1
             do while(p>0)
                child=grid(d)%patch%elts(p+1)%children(corner-4)
                if(grid(d)%patch%elts(p+1)%level<l) then
                   p=child;cycle
                end if
                if(grid(d)%patch%elts(p+1)%level>l) exit
                if(child>0) then
                   if(corner==IMINUSJPLUS.or.corner==IPLUSJMINUS) capacity=capacity+1
                   if(corner==IJMINUS) capacity=capacity+2
                end if
                p=child
             end do
          end do
       end do
    end do
    allocate(program%outer(capacity))
    do l=level_start-1,level_end-1

       do d=1,size(grid)
          do pi=1,grid(d)%lev(l)%length
             p=grid(d)%lev(l)%elts(pi)
             if (.not. any(grid(d)%patch%elts(p+1)%children>0)) cycle
             call get_offs_Domain(grid(d),p,op,dp0)
             do c=1,N_CHDRN
                child=grid(d)%patch%elts(p+1)%children(c)
                if (child==0) cycle
                call get_offs_Domain(grid(d),child,oc,dc)
                do j=1,PATCH_SIZE/2+1
                   jc=2*(j-1)
                   jp=j-1+chd_offs(2,c)
                   do i=1,PATCH_SIZE/2+1
                      ic=2*(i-1)
                      ip=i-1+chd_offs(1,c)
                      do e=RT,UP
                         item=Transform_Outer(level=l)
                         item%component=e+1
                         if(i<=PATCH_SIZE/2.and.j<=PATCH_SIZE/2) item%cell=idx(ic,jc,oc,dc)+1
                         a=idx(ic+end_pt(1,1,e+1),jc+end_pt(2,1,e+1),oc,dc)
                         b=idx(ic+end_pt(1,2,e+1),jc+end_pt(2,2,e+1),oc,dc)
                         item%target=[edge_slot(d,EDGE*idx(ip,jp,op,dp0)+e), &
                              edge_slot(d,EDGE*a+e),edge_slot(d,EDGE*b+e)]
                         item%weight=base_weight+[(grid(d)%I_u_wgt%elts(b+1)%enc(q),q=1,9)]
                         item%source(1)=item%target(1)
                         item%source(2)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+3)
                         item%source(3)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+4)
                         item%source(4)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+6)
                         item%source(5)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+1)
                         item%source(6)=side(ip+opp_no(1,1,e+1),jp+opp_no(2,1,e+1),hex_s_offs(e+1)+2)
                         item%source(7)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+3)
                         item%source(8)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+4)
                         item%source(9)=side(ip+opp_no(1,1,e+1),jp+opp_no(2,1,e+1),hex_s_offs(e+1)+5)
                         item%source(10)=side(ip+opp_no(1,2,e+1),jp+opp_no(2,2,e+1),hex_s_offs(e+1)+5)
                         item%source(11)=side(ip+end_pt(1,2,e+1),jp+end_pt(2,2,e+1),hex_s_offs(e+1)+6)
                         item%source(12)=side(ip+end_pt(1,1,e+1),jp+end_pt(2,1,e+1),hex_s_offs(e+1)+1)
                         item%source(13)=side(ip+opp_no(1,2,e+1),jp+opp_no(2,2,e+1),hex_s_offs(e+1)+2)
                         call append_operation
                      end do
                   end do
                end do
             end do
          end do
          ! Preserve the original per-Domain correction order, after all
          ! regular plus-side writes on that Domain and before the next one.
          do corner=NORTHEAST,NORTHWEST
             if (.not. grid(d)%penta(corner)) cycle
             p=1
             do while (p>0)
                child=grid(d)%patch%elts(p+1)%children(corner-4)
                if (grid(d)%patch%elts(p+1)%level<l) then
                   p=child
                   cycle
                end if
                if (grid(d)%patch%elts(p+1)%level>l) exit
                if (child>0) then
                   call get_offs_Domain(grid(d),p,op,dp0)
                   call get_offs_Domain(grid(d),child,oc,dc)
                   if (corner==IMINUSJPLUS) then
                      item=Transform_Outer(level=l,kind=1)
                      a=idx(0,LAST-1,oc,dc)
                      b=idx(0,LAST,oc,dc)
                      item%target(2:3)=[edge_slot(d,EDGE*a+UP),edge_slot(d,EDGE*b+UP)]
                      item%source(1:2)=[edge_slot(d,EDGE*idx(0,PATCH_SIZE,op,dp0)+UP), &
                           edge_slot(d,EDGE*idx(-1,PATCH_SIZE,op,dp0)+RT)]
                      item%weight(1)=base_weight(8)+grid(d)%I_u_wgt%elts(b+1)%enc(8)
                      call append_operation
                   else if (corner==IPLUSJMINUS) then
                      item=Transform_Outer(level=l,kind=1)
                      a=idx(LAST-1,0,oc,dc)
                      b=idx(LAST,0,oc,dc)
                      item%target(2:3)=[edge_slot(d,EDGE*a+RT),edge_slot(d,EDGE*b+RT)]
                      item%source(1:2)=[edge_slot(d,EDGE*idx(PATCH_SIZE,0,op,dp0)+RT), &
                           edge_slot(d,EDGE*idx(PATCH_SIZE,-1,op,dp0)+UP)]
                      item%weight(1)=-(base_weight(7)+grid(d)%I_u_wgt%elts(b+1)%enc(7))
                      call append_operation
                   else if (corner==IJMINUS) then
                      do mode=2,3
                         item=Transform_Outer(level=l,kind=mode)
                         e=merge(UP,RT,mode==2)
                         a=idx(0,0,oc,dc)
                         b=idx(merge(0,1,mode==2),merge(1,0,mode==2),oc,dc)
                         item%target(2:3)=[edge_slot(d,EDGE*a+e),edge_slot(d,EDGE*b+e)]
                         t=idx(end_pt(1,2,e+1),end_pt(2,2,e+1),oc,dc)
                         q=merge(9,6,mode==2)
                         item%weight(1)=base_weight(q)+grid(d)%I_u_wgt%elts(t+1)%enc(q)
                         if (mode==2) then
                            item%source(1:2)=[edge_slot(d,EDGE*idx(0,-1,op,dp0)+UP), &
                                 edge_slot(d,EDGE*idx(-1,-1,op,dp0)+RT)]
                            item%source(3)=side(end_pt(1,1,UP+1),end_pt(2,1,UP+1),hex_s_offs(UP+1)+1)
                            item%source(4)=side(opp_no(1,2,UP+1),opp_no(2,2,UP+1),hex_s_offs(UP+1)+2)
                         else
                            item%source(1:2)=[edge_slot(d,EDGE*idx(-1,-1,op,dp0)+RT), &
                                 edge_slot(d,EDGE*idx(-1,0,op,dp0)+RT)]
                            item%source(3)=side(opp_no(1,1,RT+1),opp_no(2,1,RT+1),hex_s_offs(RT+1)+2)
                            item%source(4)=side(end_pt(1,1,RT+1),end_pt(2,1,RT+1),hex_s_offs(RT+1)+3)
                         end if
                         call append_operation
                      end do
                   end if
                end if
                p=child
             end do
          end do
       end do

    end do
    if(noperation/=capacity) call die('outer tape count')
  contains
    integer function side(i,j,s) result(address)
      integer, intent(in) :: i,j,s
      address=edge_slot(d,ed_idx(i,j,hex_sides(:,s),op,dp0))
    end function
    subroutine append_operation
      noperation=noperation+1
      if (noperation>size(program%outer)) call die('outer tape capacity')
      item%global_domain=glo_id(rank+1,d)+1
      if(item%kind==0) then
         item%forward_active=merge(1,0,i<=PATCH_SIZE/2.and.j<=PATCH_SIZE/2.and. &
              grid(d)%mask_e%elts(abs(item%target(3)))>=ADJZONE)
      else
         item%forward_active=merge(1,0,grid(d)%mask_e%elts(abs(item%target(2)))>=ADJZONE)
      end if
      program%outer(noperation)=item
    end subroutine
  end subroutine


  subroutine check_forward_transform(input,expected,base,ncolumns,first_column,first_level)
    real(dp),intent(in) :: input(:),expected(:)
    integer,intent(in) :: base(:),ncolumns,first_column,first_level
    type(Transform_Geometry) :: program
    real(dp),allocatable :: solution(:,:),wavelet(:,:)
    integer :: d,n,l,i,k,j,a,t
    integer,allocatable :: inner_index(:)
    character(8) :: setting
    call get_environment_variable('WAVETRISK_TRANSFORM_KERNEL_CHECK',setting)
    if(trim(setting)/='1'.or.first_level/=level_start) return
    call compile_transform_geometry(program)
    do d=1,size(grid)
       n=(base(d+1)-base(d))/ncolumns
       allocate(solution(n,ncolumns),wavelet(n,ncolumns),inner_index(n/EDGE));inner_index=0
       do i=1,size(program%inner)
          if(program%inner(i)%d/=d) cycle
          inner_index(program%inner(i)%cell)=i
       end do
       solution=reshape(input(base(d)+1:base(d+1)),[n,ncolumns]);wavelet=0.0_dp
       do l=level_end-1,first_level,-1
          do i=1,size(program%outer)
             associate(op=>program%outer(i))
             if(op%level/=l.or.op%global_domain/=glo_id(rank+1,d)+1) cycle
             call execute_transform_outer(op,solution,wavelet,.false.)
             if(op%kind==0.and.op%cell>0.and.op%component==EDGE) then
                t=inner_index(op%cell)
                if(t<1) call die('forward inner cell missing')
                call execute_transform_inner(program%inner(t),solution,wavelet,.false.)
             end if
             end associate
          end do
       end do
       do k=1,ncolumns
          if(first_column+k-1<1.or.first_column+k-1>zlevels) cycle
          do j=1,n
             a=base(d)+(k-1)*n+j
             if(transfer(wavelet(j,k),0_int64)/=transfer(expected(a),0_int64)) then
                write(*,*) 'forward transform difference ',rank,d,j,k,wavelet(j,k),expected(a)
                call die('forward transform kernel differs')
             end if
          end do
       end do
       deallocate(solution,wavelet,inner_index)
    end do
  end subroutine

  subroutine execute_transform_scalar(op,solution,wavelet,inverse)
    type(Transform_Scalar),intent(in) :: op
    real(dp),intent(inout) :: solution(:,:),wavelet(:,:)
    logical,intent(in) :: inverse
    real(dp) :: s,x
    integer :: k,q
    if(inverse) then
       if(op%inverse_active==0) return
    else
       if(op%phase/=2.or.op%forward_active==0) return
    end if
    do k=1,size(solution,2)
       s=0.0_dp
       do q=1,op%n
          if(inverse.and.op%phase==1) then
             s=s+wavelet(op%source(q),k)*op%weight(q)
          else
             s=s+solution(op%source(q),k)*op%weight(q)
          end if
       end do
       if(inverse) then
          if(op%phase==1) then
             x=solution(op%base,k)-s*op%inverse_area
          else
             x=s*op%inverse_area+wavelet(op%base,k)
          end if
          solution(op%target,k)=x
       else
          wavelet(op%target,k)=solution(op%target,k)-s*op%inverse_area
       end if
    end do
  end subroutine

  subroutine execute_transform_inner(op,solution,wavelet,inverse)
    type(Transform_Inner),intent(in) :: op
    real(dp),intent(inout) :: solution(:,:),wavelet(:,:)
    logical,intent(in) :: inverse
    integer,parameter :: inputs(18)=[1,2,3,5,6,7,10,11,14,15,18,19,22,23,26,27,30,31]
    real(dp) :: u(32),g(32),curl(LORT:UPLT),x(6)
    integer :: k,q,t,s
    g=op%geometry
    do k=1,size(solution,2)
       do q=1,size(inputs)
          t=inputs(q);u(t)=solution(op%source(t),k)
       end do
       do t=LORT,UPLT
          s=4*(t-LORT)
          curl(t)=(u(s+1)*g(s+1)+u(s+2)*g(s+2)+u(s+3)*g(s+3))/g(s+4)
       end do
       x(1)=(g(9)*curl(LORT)-u(10)*g(10)-u(11)*g(11))/g(12)
       x(2)=(g(13)*curl(LORT)-u(14)*g(14)-u(15)*g(15))/g(16)
       x(3)=(g(17)*curl(LORT)-u(18)*g(18)-u(19)*g(19))/g(20)
       x(4)=(g(21)*curl(UPLT)-u(22)*g(22)-u(23)*g(23))/g(24)
       x(5)=(g(25)*curl(UPLT)-u(26)*g(26)-u(27)*g(27))/g(28)
       x(6)=(g(29)*curl(UPLT)-u(30)*g(30)-u(31)*g(31))/g(32)
       do q=1,6
          if(inverse) then
             solution(op%target(q),k)=x(q)+wavelet(op%target(q),k)
          else if(op%forward_active(q)/=0) then
             wavelet(op%target(q),k)=solution(op%target(q),k)-x(q)
          end if
       end do
    end do
  end subroutine

  subroutine execute_transform_outer(op,solution,wavelet,inverse)
    type(Transform_Outer),intent(in) :: op
    real(dp),intent(inout) :: solution(:,:),wavelet(:,:)
    logical,intent(in) :: inverse
    real(dp) :: s(13),x(9),correction,first,second
    integer :: k,q,a,b
    if(.not.inverse.and.op%forward_active==0) return
    a=op%target(2);b=op%target(3)
    do k=1,size(solution,2)
       do q=1,13
          if(op%source(q)==0) exit
          s(q)=solution(op%source(q),k)
       end do
       if(op%kind==0) then
          x(1:5)=s(1:5)
          x(6:9)=[s(6)-s(7),s(8)-s(9),s(10)-s(11),s(12)-s(13)]
          if(inverse) then
             second=sum(op%weight*x)+wavelet(b,k)
             solution(b,k)=second
             first=2.0_dp*solution(op%target(1),k)-second
             solution(a,k)=first
          else
             wavelet(b,k)=solution(b,k)-sum(op%weight*x)
             wavelet(a,k)=-wavelet(b,k)
          end if
       else
          select case(op%kind)
          case(1)
             correction=op%weight(1)*(s(1)+s(2))
          case(2)
             correction=op%weight(1)*((-s(1)-(-s(2)))-(s(3)-s(4)))
          case(3)
             correction=op%weight(1)*(s(1)+s(2)-(s(3)-s(4)))
          case default
             error stop 'transform outer operation kind'
          end select
          if(inverse) then
             solution(a,k)=solution(a,k)-correction
             solution(b,k)=solution(b,k)+correction
          else
             wavelet(a,k)=wavelet(a,k)+correction
             wavelet(b,k)=wavelet(b,k)-correction
          end if
       end if
    end do
  end subroutine
end module
