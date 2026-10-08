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
 'src/comm_mpi.f90': [],
}
BLOCK = {
 'src/time_integr.f90': [('prepare_block_multistage_tendency','TENDENCY')],
 'src/multi_level.f90': [('block_tendency_compatibility_ml','COMPAT'),('prepare_native_velocity_programs','VELOCITY_PLAN'),('prepare_native_mass','MASS_PLAN')],
 'src/parallel_block_inverse.f90': [('prepare_native_inverse','INVERSE_SETUP'),('publish_canonical_inverse203','PUBLICATION'),('install_inverse203','INSTALL'),('finish_exchange_values','INVERSE_MPI')],
 'src/parallel_block_mpi.f90': [('activate_block_native_inverse_transform','NATIVE_INVERSE'),('refresh_parallel_block_candidate_boundary_state','PROVISIONAL'),('evaluate_candidate_block_velocity_recomposition','RECOMPOSE'),('build_parallel_block_catalog','CATALOG'),('migrate_blocks','MIGRATE'),('retain_post_grid_change_block_reconstruction','POSTGRID'),('refresh_block_sol_ghosts','GHOST'),('refresh_block_sol_wav_coeff_ghosts','GHOST'),('prepare_block_native_wavelet_compression','COMPRESSION'),('activate_block_native_wavelet_compression','COMPRESSION'),('validate_candidate_block_scalar_wavelets','FORWARD_SCALAR'),('validate_candidate_block_outer_vector_wavelets','FORWARD_VECTOR'),('refresh_parallel_block_domain_prognostic_state','ACCEPTED')],
 'src/parallel_block_build.f90': [('build_source_blocks','SOURCE')],
 'src/parallel_block.f90': [('ensure_local_block_hydrostatic_state','HYDRO'),('install_local_blocks','STORE_INSTALL'),('retain_local_block_tendency_in_accumulator','ACCUM_RETAIN')],
}


BLOCK['src/parallel_block_mpi.f90'] += [
 ('begin_block_scalar_divergence_capture','CAPTURE_BEGIN'),
 ('finalize_block_scalar_divergence_capture','CAPTURE_PUBLISH'),
 ('prepare_block_native_velocity_remainder','VELOCITY_REMAINDER'),
 ('begin_block_domain_multistage_candidate_stage','RK_ASSEMBLE'),
 ('evaluate_candidate_block_scalar_restriction','SCALAR_RESTRICT'),
 ('exchange_block_scalar_restriction_ghosts','RESTRICT_GHOST'),
 ('exchange_temperature_edges','TEMP_EDGES'),
 ('evaluate_candidate_block_exner_shadow','EXNER_TRANSACTION'),
 ('refresh_candidate_block_edge_metrics','EDGE_METRIC'),
 ('exchange_domain_to_block_payloads','DOMAIN_PAYLOAD'),
 ('build_block_migration_manifest','MIG_MANIFEST'),
 ('check_block_migration_manifest','MIG_MANIFEST'),
 ('exchange_block_migration_payloads','MIG_PAYLOAD'),
 ('exchange_block_migration_sizes','MIG_SIZES'),
 ('build_block_ghost_exchange_plan','GHOST_PLAN'),
 ('build_block_writeback_plan','WRITEBACK_PLAN'),
 ('prepare_parallel_block_production_state','PRODUCTION_PREP'),
]


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



def scoped_inject(text,routine,old,new):
    start,end=bounds(text,routine)
    return text[:start]+inject(text[start:end],old,new)+text[end:]


def mirror_scopes(text,routine,function_prefix,mapping):
    """Mirror existing batch/level transaction boundaries, never sample loops."""
    start,end=bounds(text,routine);body=text[start:end]
    counts={key:[0,0] for key in mapping}
    pattern=r'^([ \t]*)call\s+'+function_prefix+r'_(enter|leave)\s*\((.*?)\)[ \t]*$'
    def add(m):
        phase=re.match(r'[ \t\n&]*([A-Z][A-Z_0-9]+)',m.group(3))
        if not phase or phase.group(1) not in mapping:return m.group(0)
        key=phase.group(1);side=m.group(2);counts[key][0 if side=='enter' else 1]+=1
        call=m.group(1)+'call cost_'+side+'(C_'+mapping[key]+')'
        return call+'\n'+m.group(0) if side=='enter' else m.group(0)+'\n'+call
    body=re.sub(pattern,add,body,flags=re.M|re.S)
    if any(a==0 or a!=b for a,b in counts.values()):raise ValueError('Batch scope anchors differ: '+routine+str(counts))
    return text[:start]+body+text[end:]


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
        path=root/'src/parallel_block_mpi.f90';text=path.read_text()
        text=mirror_scopes(text,'evaluate_candidate_block_scalar_restriction','block_profile',{'BLOCK_PROFILE_RESTRICTION_KERNEL':'SCALAR_LOCAL'})
        text=mirror_scopes(text,'exchange_block_scalar_restriction_ghosts','block_profile',{'BLOCK_PROFILE_RESTRICTION_PACK':'RESTRICT_PACK','BLOCK_PROFILE_RESTRICTION_MPI':'RESTRICT_MPI','BLOCK_PROFILE_RESTRICTION_INSTALL':'RESTRICT_INSTALL'})
        text=mirror_scopes(text,'finalize_block_scalar_divergence_capture','detail',{'DP_SCALAR_MPI':'CAPTURE_MPI'})
        text=mirror_scopes(text,'exchange_temperature_edges','detail',{'DP_SCALAR_MPI':'TEMP_WAIT'})
        text=mirror_scopes(text,'migrate_blocks','detail',{'DP_PACK':'MIG_PACK'})
        text=scoped_inject(text,'migrate_blocks','    do i = 1, manifest%n_recv\n','    call cost_enter(C_MIG_UNPACK)\n    do i = 1, manifest%n_recv\n')
        text=scoped_inject(text,'migrate_blocks','    if (pos /= n_recv_byte) then\n','    call cost_leave(C_MIG_UNPACK)\n    if (pos /= n_recv_byte) then\n')
        text=scoped_inject(text,'evaluate_candidate_block_exner_shadow',
              '    call apply_local_block_tendency_kernel( &\n         compute_block_exner_difference_kernel,kernel)\n',
              '    call cost_enter(C_EXNER_KERNEL)\n    call apply_local_block_tendency_kernel( &\n         compute_block_exner_difference_kernel,kernel)\n    call cost_leave(C_EXNER_KERNEL)\n')
        text=scoped_inject(text,'exchange_block_migration_payloads','    n_send_byte = int(manifest%total_send_nbyte)\n',
              '    call cost_work(C_MIG_PAYLOAD,manifest%total_send_nbyte,manifest%total_recv_nbyte,1_int64)\n    n_send_byte = int(manifest%total_send_nbyte)\n')
        text=scoped_inject(text,'build_block_ghost_exchange_plan','       ghost_exchange_plan%reuse_count = &\n',
              '       call cost_work(C_GHOST_PLAN,0_int64,1_int64,0_int64)\n       ghost_exchange_plan%reuse_count = &\n')
        text=scoped_inject(text,'build_block_ghost_exchange_plan','    call clear_block_ghost_exchange_plan\n',
              '    call cost_work(C_GHOST_PLAN,1_int64,0_int64,0_int64)\n    call clear_block_ghost_exchange_plan\n')
        text=scoped_inject(text,'build_block_writeback_plan','       block_writeback_plan_reuse_count = &\n',
              '       call cost_work(C_WRITEBACK_PLAN,0_int64,1_int64,0_int64)\n       block_writeback_plan_reuse_count = &\n')
        text=scoped_inject(text,'build_block_writeback_plan','    call clear_block_writeback_plan\n',
              '    call cost_work(C_WRITEBACK_PLAN,1_int64,0_int64,0_int64)\n    call clear_block_writeback_plan\n')
        text=scoped_inject(text,'refresh_candidate_block_edge_metrics','    if (storage_current) then\n       call cost_leave(C_EDGE_METRIC)\n       return\n    end if\n',
              '    if (storage_current) then\n       call cost_work(C_EDGE_METRIC,0_int64,1_int64,0_int64)\n       call cost_leave(C_EDGE_METRIC)\n       return\n    end if\n    call cost_work(C_EDGE_METRIC,1_int64,0_int64,0_int64)\n')
        text=scoped_inject(text,'exchange_block_migration_payloads','    call MPI_Alltoallv( &\n',
              '    call cost_enter(C_MIG_TRANSFER_MPI)\n    call MPI_Alltoallv( &\n')
        text=scoped_inject(text,'exchange_block_migration_payloads','    call check_mpi(ierr,"MPI_Alltoallv packed block payloads")\n',
              '    call check_mpi(ierr,"MPI_Alltoallv packed block payloads")\n    call cost_leave(C_MIG_TRANSFER_MPI)\n    call cost_enter(C_MIG_CHECKSUM)\n')
        text=scoped_inject(text,'exchange_block_migration_payloads',
              '    call MPI_Allreduce( &\n         local_send_checksum, global_send_checksum, 1, &\n',
              '    call cost_leave(C_MIG_CHECKSUM)\n    call cost_enter(C_MIG_PROOF_MPI)\n    call MPI_Allreduce( &\n         local_send_checksum, global_send_checksum, 1, &\n')
        text=scoped_inject(text,'exchange_block_migration_payloads','    call check_mpi(ierr,"MPI_Allreduce transferred payload bytes")\n',
              '    call check_mpi(ierr,"MPI_Allreduce transferred payload bytes")\n    call cost_leave(C_MIG_PROOF_MPI)\n')
        text=scoped_inject(text,'exchange_block_scalar_restriction_ghosts','    call complete_sparse_payload(request_count)\n',
              '    call complete_sparse_payload(request_count)\n    call cost_work(C_RESTRICT_MPI,int(request_count,int64), &\n         sum(int(merge(block_scalar_restriction_exchange%ghost_send_count, &\n         block_scalar_restriction_exchange%ghost_dynamic_send_count,full_payload),int64)) - &\n         int(merge(block_scalar_restriction_exchange%ghost_send_count(rank+1), &\n         block_scalar_restriction_exchange%ghost_dynamic_send_count(rank+1),full_payload),int64), &\n         sum(int(merge(block_scalar_restriction_exchange%ghost_recv_count, &\n         block_scalar_restriction_exchange%ghost_dynamic_recv_count,full_payload),int64)) - &\n         int(merge(block_scalar_restriction_exchange%ghost_recv_count(rank+1), &\n         block_scalar_restriction_exchange%ghost_dynamic_recv_count(rank+1),full_payload),int64))\n')
        text=scoped_inject(text,'exchange_temperature_edges',
              '      call record_parallel_block_profile_volume(BLOCK_PROFILE_TEMPERATURE_EDGES, &\n',
              '      call cost_work(C_TEMP_EDGES,int(nrequest,int64), &\n           int(zlevels,int64)*(sum(int(plan%service_count,int64))-int(plan%service_count(rank+1),int64)), &\n           int(zlevels,int64)*(sum(int(plan%request_count,int64))-int(plan%request_count(rank+1),int64)))\n      call record_parallel_block_profile_volume(BLOCK_PROFILE_TEMPERATURE_EDGES, &\n')
        path.write_text(text)
    return sorted(changed)
