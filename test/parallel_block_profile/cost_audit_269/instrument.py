"""Add bounded observers to diagnostic copies; exact controls remain untouched."""
from pathlib import Path
import re
import shutil

HERE = Path(__file__).resolve().parent
COMMON = {
 'src/main.f90': [('time_step','STEP')],
 'src/time_integr.f90': [('RK4','RK')],
 'src/adapt.f90': [('WT_after_step','WT'),('adapt','ADAPT')],
 'src/wavelet.f90': [('inverse_wavelet_transform','DOMAIN_INVERSE')],
 'src/multi_level.f90': [],
 'src/diagnostics.f90': [('cal_surf_press','PRESSURE')],
 'src/comm_mpi.f90': [(name,'BOUNDARY') for name in ('update_bdry_0','update_bdry_1','update_bdry_2','update_bdry1_0','update_bdry1_1','update_bdry1_2')],
}
BLOCK = {
 'src/time_integr.f90': [('prepare_block_multistage_tendency','TENDENCY')],
 'src/multi_level.f90': [('compute_native_velocity_source','VSOURCE'),('compute_native_velocity_gradient','VGRADIENT')],
 'src/parallel_block_inverse.f90': [('prepare_native_inverse','INVERSE_SETUP'),('publish_canonical_inverse203','PUBLICATION'),('install_inverse203','INSTALL'),('finish_exchange_values','INVERSE_MPI')],
 'src/parallel_block_mpi.f90': [('activate_block_native_inverse_transform','NATIVE_INVERSE'),('refresh_parallel_block_candidate_boundary_state','PROVISIONAL'),('evaluate_candidate_block_velocity_recomposition','RECOMPOSE'),('build_parallel_block_catalog','CATALOG'),('migrate_blocks','MIGRATE'),('retain_post_grid_change_block_reconstruction','POSTGRID'),('refresh_block_sol_ghosts','GHOST'),('refresh_block_sol_wav_coeff_ghosts','GHOST'),('prepare_block_native_wavelet_compression','COMPRESSION'),('activate_block_native_wavelet_compression','COMPRESSION'),('validate_candidate_block_scalar_wavelets','FORWARD_SCALAR'),('validate_candidate_block_outer_vector_wavelets','FORWARD_VECTOR'),('refresh_parallel_block_domain_prognostic_state','ACCEPTED')],
 'src/parallel_block_build.f90': [('build_source_blocks','SOURCE')],
 'src/parallel_block.f90': [('ensure_local_block_hydrostatic_state','HYDRO')],
}

def bounds(text, name):
    starts=list(re.finditer(r'^[ \t]*subroutine\s+'+re.escape(name)+r'\b[^\n]*',text,re.M|re.I))
    if len(starts)!=1: raise ValueError('Require one subroutine: '+name)
    start=starts[0].start();depth=0
    for m in re.finditer(r'^[ \t]*(subroutine\s+\w+|end\s+subroutine\b)[^\n]*',text[start:],re.M|re.I):
        depth += -1 if m.group(1).lower().startswith('end') else 1
        if depth==0:return start,start+m.end()
    raise ValueError('Unclosed subroutine: '+name)


def wrap(text,name,phase):
    start,end=bounds(text,name)
    body=text[start:end]
    cut=re.search(r'^[ \t]*contains\s*$',body,re.M|re.I)
    tail=re.search(r'^[ \t]*end\s+subroutine\b',body,re.M|re.I) if cut is None else cut
    if tail is None:raise ValueError('Missing routine tail: '+name)
    outer,rest=body[:tail.start()],body[tail.start():]
    lines=outer.splitlines(keepends=True)
    executable=None
    for i,line in enumerate(lines):
        stripped=line.strip()
        if not stripped or stripped.startswith(('!','#')) or '::' in stripped:continue
        if re.match(r'(call\b|if\s*\(|do\b|select\b|where\s*\(|allocate\s*\(|deallocate\s*\()',stripped,re.I) or re.match(r'\w+(?:%\w+)*(?:\([^=]*\))?\s*=(?!=)',stripped):
            executable=i;break
    if executable is None:raise ValueError('Cannot locate first executable: '+name)
    enter='cost_begin_step' if phase=='STEP' else 'cost_enter(C_'+phase+')'
    leave='cost_end_step' if phase=='STEP' else 'cost_leave(C_'+phase+')'
    lines.insert(executable,'    call '+enter+'\n')
    outer=''.join(lines)
    def close_return(m):
        indent,condition=m.group(1),m.group(2)
        if condition:return indent+condition+' then\n'+indent+'   call '+leave+'\n'+indent+'   return\n'+indent+'end if'
        return indent+'call '+leave+'\n'+indent+'return'
    outer=re.sub(r'^([ \t]*)(if\s*\([^\n]*\))?\s*return\s*$',close_return,outer,flags=re.M|re.I)
    outer+='    call '+leave+'\n'
    return text[:start]+outer+rest+text[end:]


def inject(text,old,new):
    if text.count(old)!=1:raise ValueError('Instrumentation anchor differs: '+old)
    return text.replace(old,new)


def instrument(root,block):
    root=Path(root)
    specs={name:list(rows) for name,rows in COMMON.items()}
    specs['src/multi_level.f90']=[] if block else [('trend_ml','TENDENCY')]
    if block:
        for name,rows in BLOCK.items():specs.setdefault(name,[]).extend(rows)
    changed=[]
    for name,rows in specs.items():
        if not rows:continue
        path=root/name;text=path.read_text()
        module=re.search(r'^module\s+\w+[^\n]*',text,re.M|re.I)
        if module is None:raise ValueError('Module declaration missing: '+name)
        text=text[:module.end()]+'\n  use cost_observer_mod'+text[module.end():]
        for symbol,phase in rows:text=wrap(text,symbol,phase)
        path.write_text(text);changed.append(name)
    path=root/'test/climate/climate.f90'
    text=path.read_text()
    if not text.startswith('program climate\n'):raise ValueError('Climate program anchor differs')
    text=text.replace('program climate\n','program climate\n  use cost_observer_mod, only : cost_flush\n',1)
    text=inject(text,'  call finalize\n','  call cost_flush(rank)\n  call finalize\n')
    path.write_text(text);changed.append('test/climate/climate.f90')
    shutil.copy2(HERE/'cost_observer.f90',root/'src/cost_observer.f90');changed.append('src/cost_observer.f90')
    make=root/'Makefile'
    text=inject(make.read_text(),'SRC = kind.f90 \\\n','SRC = cost_observer.f90 \\\n      kind.f90 \\\n')
    objects=' '.join('$(BUILD_DIR)/'+Path(name).stem+'.o' for name in changed if name.startswith('src/') and name!='src/cost_observer.f90')
    text+='\n'+objects+' $(BUILD_DIR)/test.o: $(BUILD_DIR)/cost_observer.o\n'
    make.write_text(text);changed.append('Makefile')
    if block:
        path=root/'src/parallel_block.f90';text=path.read_text()
        text=inject(text,'  block_hydrostatic(local_index)%ready = .true.\n',
                    '  call cost_work(C_HYDRO,int(n_surface,int64),int(n_column,int64),1_int64)\n  block_hydrostatic(local_index)%ready = .true.\n')
        path.write_text(text)
        path=root/'src/parallel_block_inverse.f90';text=path.read_text()
        text=inject(text,'    if (plan_generation /= generation) then\n','    if (plan_generation /= generation) then\n       call cost_work(C_INVERSE_SETUP,1_int64,0_int64,0_int64)\n')
        text=inject(text,'    installation203_count=installation203_count+1_int64\n',
                    '    call cost_work(C_INSTALL,int(size(interior_plan%slot),int64), &\n         int(size(boundary_plan%slot),int64),1_int64)\n    installation203_count=installation203_count+1_int64\n')
        text=inject(text,'    scaling%bdry_uptodate=.false.\n',
                    '    call cost_work(C_PUBLICATION,int(size(publication_patch)*PATCH_SIZE**2,int64), &\n         int(size(publication_patch),int64)*int(size(publication_tile),int64),1_int64)\n    scaling%bdry_uptodate=.false.\n')
        text=inject(text,"    call detail_leave(DP_INVERSE_POST)\n",
                    "    call cost_work(C_INVERSE_MPI,int(nrequest,int64), &\n         int(n,int64)*sum(int(sc,int64))-int(n,int64)*int(sc(rank+1),int64), &\n         int(n,int64)*sum(int(rc,int64))-int(n,int64)*int(rc(rank+1),int64))\n    call detail_leave(DP_INVERSE_POST)\n")
        path.write_text(text)
    return sorted(changed)
