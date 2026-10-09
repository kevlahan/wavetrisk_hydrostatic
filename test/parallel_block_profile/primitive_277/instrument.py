"""Integrate fused primitive arithmetic without adding an intermediate MPI route."""
from pathlib import Path

def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Source anchor differs: '+old[:90])
    return text.replace(old,new,1)

def instrument(root):
    path=Path(root)/'src/multi_level.f90';text=path.read_text()
    text=once(text,'       post_step1, Qperp, scalar_trend, scalar_trend_pair, step1',
              '       post_step1, Qperp, scalar_trend, scalar_trend_pair, step1, prepare_primitive277, finish_primitive277')
    text=once(text,'    call prepare_native_mass\n','    call prepare_native_mass\n    call prepare_primitive277(native_velocity_plan_generation())\n')
    text=once(text,'    native_velocity_transaction=.false.','    call finish_primitive277\n    native_velocity_transaction=.false.')
    path.write_text(text)
    path=Path(root)/'src/ops.f90';text=path.read_text()
    text=once(text,'module ops_mod\n\n','module ops_mod\n  use, intrinsic :: iso_fortran_env, only : int64\n  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite\n  use arch_mod, only : rank,n_process\n  use shared_mod, only : N_GLO_DOMAIN\n')
    text=once(text,'  public :: comp_offs3','  public :: comp_offs3, prepare_primitive277, finish_primitive277')
    text=once(text,'\ncontains\n',Path(__file__).with_name('declarations.f90').read_text()+'\ncontains\n')
    text=once(text,'    logical :: S_bdry, W_bdry, skip_temperature',
              '    logical :: S_bdry, W_bdry, skip_temperature, fused_patch\n    real(dp) :: primitive_patch(19,PATCH_SIZE**2), primitive_x(16)\n    integer :: primitive_nn,primitive_ids(6),primitive_i,primitive_j')
    anchor='    id = offs(0)\n'
    setup='''    fused_patch=itype==0.and.mass_transaction
    if (fused_patch) then
       if (.not.allocated(primitive277)) error stop 'Stage277 cache missing'
       if (.not.allocated(primitive277(dom%id+1)%patch(p+1)%geometry)) error stop 'Stage277 patch cache missing'
       if (p277_check) p277_layer_nodes(zlev)=p277_layer_nodes(zlev)+PATCH_SIZE**2
       do primitive_j=0,PATCH_SIZE-1
          do primitive_i=0,PATCH_SIZE-1
             primitive_nn=primitive_j*PATCH_SIZE+primitive_i+1
             primitive_ids=primitive277(dom%id+1)%patch(p+1)%mass(:,primitive_nn)
             primitive_x(1:6)=q(S_MASS,zlev)%data(dom%id+1)%elts(primitive_ids)+ &
                  sol_mean(S_MASS,zlev)%data(dom%id+1)%elts(primitive_ids)
             primitive_x(7:16)=velo(primitive277(dom%id+1)%patch(p+1)%velocity(:,primitive_nn))
             call fused_primitive277(primitive_x,primitive277(dom%id+1)%patch(p+1)%geometry(:,primitive_nn), &
                  primitive_patch(:,primitive_nn))
          end do
       end do
    end if

'''
    text=once(text,anchor,setup+anchor)
    text=once(text,'      real(dp) :: pv_mass(4)','      real(dp) :: pv_mass(4)\n      integer :: pn\n      logical :: fused_node\n      real(dp) :: py(19),reference(16)')
    text=once(text,'         id_rhodz = [ id, idN, idE, idS, idW, idNE ]',
              '         fused_node=fused_patch.and.id>=offs(0).and.id<offs(0)+PATCH_SIZE**2\n         id_rhodz = [ id, idN, idE, idS, idW, idNE ]')
    anchor='         do v = scalars(1), scalars(2)\n            if (skip_temperature .and. v == S_TEMP .and. .not. any(temperature_needed)) cycle'
    first=text.index(anchor)
    text=text[:first]+text[first:].replace(anchor,anchor.replace('            if (skip_temperature','            if (fused_node.and..not.p277_check.and.v==S_MASS) cycle\n            if (skip_temperature'),1)
    start=text.index('         u_dual_RT    =',text.index('      elseif (itype == 0) then ! standard'))
    end=text.index('         ! Interpolate geopotential',start)
    old=text[start:end]
    new='''         fused_node=fused_patch.and.id>=offs(0).and.id<offs(0)+PATCH_SIZE**2
         if (fused_node) then
            pn=id-offs(0)+1
            py=primitive_patch(:,pn)
         end if
         if (.not.fused_node.or.p277_check) then
'''+old+'''         end if
         if (fused_node) then
            if (p277_check) then
               reference=[qe(EDGE*id+1:EDGE*id+3),ke(id_i), &
                    vort(TRIAG*id+1:TRIAG*id+2),vort(TRIAG*idW+LORT+1),vort(TRIAG*idS+UPLT+1), &
                    pv_W_LORT,pv_S_UPLT,u_dual_RT,u_dual_DG,u_dual_UP,u_dual_RT_W,u_dual_DG_SW,u_dual_UP_S]
               call assert_primitive277(py(4:19),reference)
            end if
            qe(EDGE*id+1:EDGE*id+3)=py(4:6)
            ke(id_i)=py(7)
            vort(TRIAG*id+1:TRIAG*id+2)=py(8:9)
            vort(TRIAG*idW+LORT+1)=py(10)
            vort(TRIAG*idS+UPLT+1)=py(11)
            pv_W_LORT=py(12);pv_S_UPLT=py(13)
            u_dual_RT=py(14);u_dual_DG=py(15);u_dual_UP=py(16)
            u_dual_RT_W=py(17);u_dual_DG_SW=py(18);u_dual_UP_S=py(19)
         end if

'''
    text=text[:start]+new+text[end:]
    start=text.index('            scalar_flux(EDGE*id+RT+1) = u_dual_RT')
    end=text.index('            if (mass_transaction.and.mass_oracle.and.v==S_MASS)',start)
    old=text[start:end]
    new='''            if (.not.fused_node.or.v/=S_MASS.or.p277_check) then
'''+old+'''            end if
            if (fused_node.and.v==S_MASS) then
               if (p277_check) call assert_primitive277(py(1:3)+physics_flux,scalar_flux(EDGE*id+1:EDGE*id+3))
               scalar_flux(EDGE*id+1:EDGE*id+3)=py(1:3)+physics_flux
            end if
'''
    text=text[:start]+new+text[end:]
    text=once(text,'end module ops_mod',Path(__file__).with_name('producer.f90').read_text()+'\nend module ops_mod')
    path.write_text(text)
    return ['src/multi_level.f90','src/ops.f90']
