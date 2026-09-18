#!/usr/bin/env python3
"""Compile generic Domain halo routes; preserve field/wire and alias order."""
from pathlib import Path
from candidate_reconstruction_209 import once
NAMES=['src/comm_mpi.f90']
DECL='''  ! Stage215 caches INTEGER routes only. Numerical buffers remain per exchange.
  type :: Boundary_Group215
     integer :: peer=0, domain=0
     type(Int_Array) :: index(AT_NODE:AT_EDGE)
  end type
  type :: Boundary_Plan215
     integer :: lo=0, hi=-1
     type(Boundary_Group215), allocatable :: send(:), receive(:)
     integer, allocatable :: sc(:,:), rc(:,:)
     type(Boundary_Plan215), pointer :: next=>null()
  end type
  type(Boundary_Plan215), pointer :: plans215=>null()
'''
HELP='''  subroutine invalidate_boundary215
    type(Boundary_Plan215), pointer :: old
    ! comm_communication_mpi completes every init/refinement route rebuild.
    ! This also handles restart reconstruction; no recycled indices survive.
    do while(associated(plans215))
       old=>plans215
       plans215=>old%next
       deallocate(old)
    end do
  end subroutine

  subroutine compile_boundary215(plan,lo,hi)
    type(Boundary_Plan215), intent(inout) :: plan
    integer, intent(in) :: lo,hi
    integer :: direction,r,ds,dd,d,g,pos,m,i,id,n,total,ng
    type(Boundary_Group215), allocatable :: groups(:)
    type(Int_Array) :: raw
    plan%lo=lo;plan%hi=hi
    allocate(plan%sc(AT_NODE:AT_EDGE,n_process),plan%rc(AT_NODE:AT_EDGE,n_process))
    plan%sc=0;plan%rc=0
    total=n_domain(rank+1)*sum(n_domain)
    do direction=1,2
       allocate(groups(total));ng=0
       do r=1,n_process
          if(r==rank+1) cycle
          do ds=1,merge(n_domain(rank+1),n_domain(r),direction==1)
             do dd=1,merge(n_domain(r),n_domain(rank+1),direction==1)
                d=merge(ds,dd,direction==1)
                g=glo_id(r,merge(dd,ds,direction==1))+1
                ng=ng+1;groups(ng)%domain=d;groups(ng)%peer=r
                do pos=AT_NODE,AT_EDGE
                   m=merge(EDGE,1,pos==AT_EDGE)
                   if(direction==1) then
                      raw=grid(d)%pack(pos,g)
                   else
                      raw=grid(d)%unpk(pos,g)
                   end if
                   n=0
                   do i=1,raw%length
                      id=abs(raw%elts(i))/m+1
                      if(grid(d)%level%elts(id)<lo.or.grid(d)%level%elts(id)>hi) cycle
                      n=n+1
                   end do
                   groups(ng)%index(pos)%length=n
                   allocate(groups(ng)%index(pos)%elts(n));n=0
                   do i=1,raw%length
                      id=abs(raw%elts(i))/m+1
                      if(grid(d)%level%elts(id)<lo.or.grid(d)%level%elts(id)>hi) cycle
                      n=n+1;groups(ng)%index(pos)%elts(n)=raw%elts(i)
                   end do
                   if(direction==1) then
                      plan%sc(pos,r)=plan%sc(pos,r)+n
                   else
                      plan%rc(pos,r)=plan%rc(pos,r)+n
                   end if
                end do
                if(all(groups(ng)%index%length==0)) then
                   do pos=AT_NODE,AT_EDGE
                      deallocate(groups(ng)%index(pos)%elts)
                   end do
                   ng=ng-1
                end if
             end do
          end do
       end do
       if(direction==1) then
          plan%send=groups(:ng)
       else
          plan%receive=groups(:ng)
       end if
       deallocate(groups)
    end do
  end subroutine

  subroutine select_boundary215(plan,lo,hi)
    type(Boundary_Plan215), pointer :: plan
    integer, intent(in) :: lo,hi
#ifdef STAGE215_AUDIT
    type(Boundary_Plan215) :: fresh
    integer :: k,pos
#endif
    plan=>plans215
    do while(associated(plan))
       if(plan%lo==lo.and.plan%hi==hi) exit
       plan=>plan%next
    end do
    if(.not.associated(plan)) then
       allocate(plan);call compile_boundary215(plan,lo,hi)
       plan%next=>plans215;plans215=>plan
    end if
#ifdef STAGE215_AUDIT
    ! Re-read the original lists on EVERY use to expose stale topology/levels.
    call compile_boundary215(fresh,lo,hi)
    if(any(fresh%sc/=plan%sc).or.any(fresh%rc/=plan%rc)) error stop 'Stage215 stale counts'
    if(size(fresh%send)/=size(plan%send).or.size(fresh%receive)/=size(plan%receive)) &
         error stop 'Stage215 stale groups'
    do k=1,size(plan%send)
       if(fresh%send(k)%domain/=plan%send(k)%domain.or.fresh%send(k)%peer/=plan%send(k)%peer) &
            error stop 'Stage215 stale send ordering'
       do pos=AT_NODE,AT_EDGE
          if(fresh%send(k)%index(pos)%length/=plan%send(k)%index(pos)%length) error stop 'Stage215 send extent'
          if(any(fresh%send(k)%index(pos)%elts/=plan%send(k)%index(pos)%elts)) error stop 'Stage215 stale send'
       end do
    end do
    do k=1,size(plan%receive)
       if(fresh%receive(k)%domain/=plan%receive(k)%domain.or.fresh%receive(k)%peer/=plan%receive(k)%peer) &
            error stop 'Stage215 stale receive ordering'
       do pos=AT_NODE,AT_EDGE
          if(fresh%receive(k)%index(pos)%length/=plan%receive(k)%index(pos)%length) error stop 'Stage215 recv extent'
          if(any(fresh%receive(k)%index(pos)%elts/=plan%receive(k)%index(pos)%elts)) error stop 'Stage215 stale recv'
       end do
    end do
#endif
  end subroutine

'''

def field_loop(rank,body):
    if rank==0:return body
    if rank==1:return '    do i1=1,size(field)\n'+body+'    end do\n'
    return '    do i2=1,size(field,2)\n       do i1=1,size(field,1)\n'+body+'       end do\n    end do\n'

def candidate(root):
    p=root/NAMES[0];s=p.read_text()
    s=s.replace('  implicit none', '''#ifdef STAGE215_AUDIT
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan,ieee_is_nan
#endif
  implicit none''',1)
    s=once(s,'\ncontains\n','\n'+DECL+'\ncontains\n'+HELP)
    s=once(s,'    call recreate_send_patch_lists()','    call recreate_send_patch_lists()\n    call invalidate_boundary215')
    # init invalidates even if an unusual caller initializes without refinement.
    s=once(s,'    call init(send_buf_i, 0)','    call invalidate_boundary215\n    call init(send_buf_i, 0)')
    for rank in range(3):
        f=['field','field(i1)','field(i1,i2)'][rank]
        for kind in ('start','finish'):
            name=f'update_bdry__{kind}1_{rank}'
            lo=s.index('  subroutine '+name);hi=s.index('  end subroutine '+name,lo)
            part=s[lo:hi]
            import re
            part=re.sub(r'^    integer ::.*\n','',part,flags=re.M)
            if rank==0:
                part=part.replace("    if (field%pos == AT_EDGE) then\n       multipl = EDGE\n    else\n       multipl = 1\n    end if\n",'')
            variables=['group215','pos','i','id']
            variables+=['d_src','r','tag','j215','count215(AT_NODE:AT_EDGE)'] if kind=='start' else ['d_dest','k']
            if rank>0:variables+=['i1']
            if rank>1:variables+=['i2']
            insert='    type(Boundary_Plan215), pointer :: plan215\n    integer :: '+', '.join(variables)+'\n'
            part=part.replace('    if (',insert+'\n    if (',1)
            if kind=='start':
                a=part.index('    send_buf%length = 0');b=part.index('    call detail_leave(DP_BDRY_ROUTE)',a)
                code='''    call select_boundary215(plan215,l_start,l_end)
    count215=0
'''+field_loop(rank,f'       count215({f}%pos)=count215({f}%pos)+1\n')+'''    do r=1,n_process
       send_lengths(r)=sum(count215*plan215%sc(:,r))
       recv_lengths(r)=sum(count215*plan215%rc(:,r))
    end do
    send_offsets(1)=0;recv_offsets(1)=0
    do r=2,n_process
       send_offsets(r)=send_offsets(r-1)+send_lengths(r-1)
       recv_offsets(r)=recv_offsets(r-1)+recv_lengths(r-1)
    end do
    send_buf%length=sum(send_lengths);recv_buf%length=sum(recv_lengths)
    if(allocated(send_buf%elts)) then
       if(size(send_buf%elts)<send_buf%length) deallocate(send_buf%elts)
    end if
    if(.not.allocated(send_buf%elts)) allocate(send_buf%elts(send_buf%length))
    j215=0
    do group215=1,size(plan215%send)
       d_src=plan215%send(group215)%domain
'''+field_loop(rank,f'''       pos={f}%pos
       do i=1,plan215%send(group215)%index(pos)%length
          id=plan215%send(group215)%index(pos)%elts(i)
          j215=j215+1;send_buf%elts(j215)={f}%data(d_src)%elts(id+1)
       end do
''')+'''    end do
    if(j215/=send_buf%length) error stop 'Stage215 pack length'
    call detail_leave(DP_BDRY_PACK)
    call detail_enter(DP_BDRY_ROUTE)
'''
                part=part[:a]+code+part[b:]
                part=once(part,'    if (recv_buf%length > 0) recv_buf%elts = 0.0_dp', '''    ! Receives exactly cover the active prefix. The tail has no consumer.
#ifdef STAGE215_AUDIT
    if(recv_buf%length>0) recv_buf%elts=ieee_value(0.0_dp,ieee_quiet_nan)
#endif''')
                part=once(part,'''    if (detail_enabled.and.recv_buf%length>0) &
         call detail_add(DC_BDRY_ZERO,int(size(recv_buf%elts),int64))''','''    ! No production receive-buffer zeroing.''')
            else:
                a=part.index('    k = 0');b=part.index('    if (k /= recv_buf%length)',a)
                code='''    call select_boundary215(plan215,l_start,l_end)
#ifdef STAGE215_AUDIT
    if(any(ieee_is_nan(recv_buf%elts(:recv_buf%length)))) error stop 'Stage215 receive hole/NaN'
#endif
    k=0
    do group215=1,size(plan215%receive)
       d_dest=plan215%receive(group215)%domain
'''+field_loop(rank,f'''       pos={f}%pos
       do i=1,plan215%receive(group215)%index(pos)%length
          id=plan215%receive(group215)%index(pos)%elts(i)
          k=k+1;{f}%data(d_dest)%elts(abs(id)+1)=recv_buf%elts(k)
          if(id<0.and.pos==AT_EDGE) {f}%data(d_dest)%elts(abs(id)+1)=-recv_buf%elts(k)
       end do
''')+'''    end do

'''
                part=part[:a]+code+part[b:]
            s=s[:lo]+part+s[hi:]
    p.write_text(s)
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('tree',type=Path);candidate(p.parse_args().tree)
