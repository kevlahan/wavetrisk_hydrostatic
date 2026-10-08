"""Apply a diagnostic-only native velocity execution-contract exporter to Stage266."""
from pathlib import Path

def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Source anchor differs: '+old[:80])
    return text.replace(old,new,1)

def instrument(root):
    path=Path(root)/'src/multi_level.f90';text=path.read_text()
    text=once(text,'module multi_level_mod\n  use parallel_block_profile_mod','module multi_level_mod\n  use arch_mod, only : audit274_catalog=>block_catalog\n  use shared_mod, only : N_GLO_DOMAIN\n  use parallel_block_profile_mod')
    # Insert immediately before module CONTAINS, not inside a routine.
    marker='\ncontains\n'
    text=once(text,marker,'\n  integer, save :: audit274_unit=-1\n  integer, allocatable, save :: audit274_samples(:)\n'+marker)
    start=text.index('  subroutine prepare_native_velocity_programs\n')
    end=text.index('  end subroutine prepare_native_velocity_programs',start)
    part=text[start:end]
    part=once(part,'    generation=native_velocity_plan_generation()\n','    call audit274_begin\n    generation=native_velocity_plan_generation()\n')
    text=text[:start]+part+text[end:]
    text=once(text,'       work%generation=generation\n       native_velocity_work(1)',
              '       work%generation=generation\n       call audit274_snapshot(d)\n       native_velocity_work(1)')
    text=once(text,'end module multi_level_mod',Path(__file__).with_name('export.f90').read_text()+'\nend module multi_level_mod')
    path.write_text(text)
    return ['src/multi_level.f90']
