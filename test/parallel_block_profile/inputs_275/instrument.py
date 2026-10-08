"""Apply a diagnostic-only native velocity execution-contract exporter to Stage266."""
from pathlib import Path

def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Source anchor differs: '+old[:80])
    return text.replace(old,new,1)

def instrument(root):
    path=Path(root)/'src/multi_level.f90';text=path.read_text()
    text=once(text,'module multi_level_mod\n  use parallel_block_profile_mod','module multi_level_mod\n  use arch_mod, only : audit275_catalog=>block_catalog\n  use shared_mod, only : N_GLO_DOMAIN\n  use parallel_block_profile_mod')
    # Insert immediately before module CONTAINS, not inside a routine.
    marker='\ncontains\n'
    text=once(text,marker,'\n  integer, save :: audit275_unit=-1\n  integer, allocatable, save :: audit275_samples(:)\n'+marker)
    start=text.index('  subroutine prepare_native_velocity_programs\n')
    end=text.index('  end subroutine prepare_native_velocity_programs',start)
    part=text[start:end]
    part=once(part,'    generation=native_velocity_plan_generation()\n','    call audit275_begin\n    generation=native_velocity_plan_generation()\n')
    text=text[:start]+part+text[end:]
    text=once(text,'       work%generation=generation\n       native_velocity_work(1)',
              '       work%generation=generation\n       call audit275_snapshot(d)\n       native_velocity_work(1)')
    text=once(text,'end module multi_level_mod',Path(__file__).with_name('export.f90').read_text()+'\nend module multi_level_mod')
    path.write_text(text)
    return ['src/multi_level.f90']

base_instrument=instrument

def instrument(root):
    changed=base_instrument(root)
    path=Path(root)/'src/multi_level.f90';text=path.read_text()
    text=once(text,'module multi_level_mod\n  use arch_mod','module multi_level_mod\n  use parallel_block_mod, only : Audit275_Block=>Block_Data, &\n       audit275_block_consumer=>apply_local_block_field_consumer\n  use arch_mod')
    text=once(text,'\ncontains\n','''
  integer, save :: audit275_probe_unit=-1,audit275_transaction=0
  integer(int64), save :: audit275_source_checks=0_int64,audit275_gradient_checks=0_int64
  real(dp), allocatable, save :: audit275_source_shadow(:,:)
  integer, parameter :: audit275_component(19)=[1,2,3,1,2,3,1,2,3,1,2,3,1,2,2,2,3,3,1]
  integer, parameter :: audit275_position(19)=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,7,7,8,9]
contains
''')
    text=once(text,'    call prepare_native_velocity_programs\n','    call audit275_transaction_begin\n    call prepare_native_velocity_programs\n')
    anchor='       call execute_velocity_sources(work%level(l),flux_view,pv_view,work%physics,work%source)'
    text=once(text,anchor,'       call audit275_source_before(d,l,flux_view,pv_view)\n'+anchor+'\n       call audit275_source_after(d)')
    start=text.index('  subroutine compute_native_velocity_gradient(q,k)')
    end=text.index('  end subroutine compute_native_velocity_gradient',start)
    part=text[start:end]
    part=once(part,'    do d=1,size(grid)','    call audit275_gradient_begin(k)\n    do d=1,size(grid)')
    part=once(part,'       work%ready(k)=.true.','       call audit275_gradient_check(d,k,q)\n       work%ready(k)=.true.')
    part=once(part,'    nullify(tendency_view)','    call audit275_gradient_end(k)\n    nullify(tendency_view)')
    text=text[:start]+part+text[end:]
    text=once(text,'end module multi_level_mod',Path(__file__).with_name('probe.f90').read_text()+'\nend module multi_level_mod')
    path.write_text(text)
    return changed
