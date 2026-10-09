"""Bounded positive-interior fused producer, isolated on accepted Stage266."""
from pathlib import Path

def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Source anchor differs: '+old[:80])
    return text.replace(old,new,1)

def instrument(root):
    path=Path(root)/'src/multi_level.f90';text=path.read_text()
    text=once(text,'module multi_level_mod\n  use parallel_block_profile_mod', '''module multi_level_mod
  use parallel_block_mod, only : Block276=>Block_Data, p276_consumer=>apply_local_block_field_consumer, &
       STORE_PATCH,STORE_BDRY,STORE_GHOST
  use arch_mod, only : block_catalog,loc_id
  use shared_mod, only : NORTH,EAST,SOUTH,WEST,NORTHEAST,NORTHWEST,SOUTHWEST,SOUTHEAST,N_GLO_DOMAIN
  use parallel_block_profile_mod''')
    text=once(text,'\ncontains\n',Path(__file__).with_name('declarations.f90').read_text()+'\ncontains\n')
    text=once(text,'    call prepare_native_velocity_programs\n','    p276_transaction=p276_transaction+1\n    call prepare_native_velocity_programs\n')
    text=once(text,'    real(dp) :: mass_start\n','    real(dp) :: mass_start\n    real(dp) :: p276_reference(9,PATCH_SIZE**2)\n')
    text=once(text,'    capture_scalar_physics = block_scalar_capture_active()','    call p276_begin(k,l)\n    capture_scalar_physics = block_scalar_capture_active()')
    old='             call step1(q=q,dom=grid(d),p=grid(d)%lev(l)%elts(j),zlev=k,itype=0, &\n                  scalar_physics=scalar_physics,native_temperature=mass_only_compatibility)'
    text=once(text,old,"             if (p276_active) then\n"+old.replace('native_temperature=mass_only_compatibility)',
              'native_temperature=mass_only_compatibility, &\n                  primitive_values=p276_reference)')+"\n             else\n"+old+"\n             end if")
    text=once(text,'             call capture_block_scalar_physics_patch( &', '''             if (p276_active) call p276_capture(q,d,grid(d)%lev(l)%elts(j),k,scalar_physics(:,:,S_MASS),p276_reference)
             call capture_block_scalar_physics_patch( &''')
    text=once(text,'    ! Retain the direct step1 flux before a coarse level', '    if (p276_active) call p276_execute(k,l)\n\n    ! Retain the direct step1 flux before a coarse level')
    text=once(text,'end module multi_level_mod',Path(__file__).with_name('producer.f90').read_text()+'\nend module multi_level_mod')
    path.write_text(text)
    path=Path(root)/'src/ops.f90';text=path.read_text()
    text=once(text,'itype, scalar_physics, native_temperature)','itype, scalar_physics, native_temperature, primitive_values)')
    text=once(text,'    logical, optional, intent(in) :: native_temperature',
              '    real(dp), optional, intent(out) :: primitive_values(9,PATCH_SIZE**2)\n    logical, optional, intent(in) :: native_temperature')
    anchor='      end if\n    end subroutine comput'
    hook="         if (present(primitive_values)) then\n            if (.not.mass_transaction) error stop 'Stage276 requires native mass'\n            if (id>=offs(0).and.id<offs(0)+PATCH_SIZE**2) then\n               primitive_values(:,id-offs(0)+1)=[native_mass(d)%flux(EDGE*id+1:EDGE*id+3), &\n                    qe(EDGE*id+1:EDGE*id+3),ke(id_i),vort(TRIAG*id+1:TRIAG*id+2)]\n            end if\n         end if\n"
    text=once(text,anchor,hook+anchor)
    path.write_text(text)
    return ['src/multi_level.f90','src/ops.f90']
