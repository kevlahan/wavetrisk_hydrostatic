"""Frozen-source, batch-only lifecycle observation; no numerical changes.

Call wrappers execute inside the original conditional. Numerical statements are
restored byte-for-byte as a preparation gate. MPI clocks operate only inside a
selected enclosing batch; no global message tape and no per-value instrumentation.
"""
import re
from pathlib import Path
GROUPS=[
 ('catalog construction',['build_parallel_block_catalog']),
 ('source block construction',['build_source_blocks']),
 ('migration transaction',['migrate_blocks']),
 ('manifest validation',['check_block_migration_manifest']),
 ('migration size routing/validation',['exchange_block_migration_sizes']),
 ('migration payload/validation',['exchange_block_migration_payloads']),
 ('payload checksum scan',[]),
 ('ghost plan construction',['build_block_ghost_exchange_plan']),
 ('writeback plan construction',['build_block_writeback_plan']),
 ('scalar capture plan',['prepare_scalar_divergence_plan']),
 ('restriction exchange plan',['prepare_block_scalar_restriction_exchange']),
 ('temperature routing plan',['prepare_temperature_edge_plan']),
 ('scalar geometry preparation',['prepare_block_scalar_capture_addresses','seed_scalar_producer_geometry']),
 ('scalar capture setup',['begin_block_scalar_divergence_capture']),
 ('scalar capture finish',['finalize_block_scalar_divergence_capture']),
 ('Domain to block transaction',['import_domain_field_family_to_blocks','exchange_domain_to_block_payloads']),
 ('block to Domain transaction',['write_block_field_family_to_domains','reconstruct_block_writeback_domain_stage']),
 ('RK tendency capture',['capture_block_domain_multistage_candidate_tendency']),
 ('RK candidate assembly',['begin_block_domain_multistage_candidate_stage']),
 ('forward vector transform',['validate_complete_block_vector_wavelets']),
 ('inverse preparation',['prepare_block_native_inverse_boundaries']),
 ('Domain RK assembly',['RK_sub_step','RK_sub_step_compatibility']),
 ('Domain basic operators',['basic_operators']),
 ('temperature boundary preparation',['compute_temperature_rk_boundary','initialize_scalar_restriction_boundary_flux']),
 ('writeback invalidation',['clear_block_writeback_plan']),
 ('selected MPI payload/metadata',['MPI_Alltoallv','MPI_Alltoall','MPI_Allgather','MPI_Allgatherv']),
 ('selected MPI reduction/dependency',['MPI_Allreduce','MPI_Reduce','MPI_Bcast','MPI_Barrier']),
 ('reserved (no per-message post clocks)',[]),
 ('selected MPI completion',['MPI_Waitall','MPI_Wait','MPI_Test','MPI_Testall'])]
REGIONS={i+1:n for i,(n,_) in enumerate(GROUPS)}
TARGETS={n.lower():i+1 for i,(_,names) in enumerate(GROUPS) for n in names}
# event = step, context, enclosing operation, RK stage, kind, object/site,
#         generation, a, b, c. Values are counts/identities, never DRAM bytes.
EVENTS={1:'plan request: object, current generation, cached generation, ready, reserved',
 2:'writeback generation installed',3:'writeback invalidation: generation, ready',
 4:'capture mode: generation, full transport, cache ready, storage rebuilt',
 5:'geometry seed request: generation',6:'import family: generation, family, vector-only, reserved',
 7:'publication family: generation, family, component, reserved',
 8:'RK tendency producer: generation, stage, stage count, reserved',
 9:'RK state update: generation, stage, stage count, reserved',
 10:'scalar cache installed: generation',11:'scalar operator consumer: generation, first level, last level, reserved',
 12:'topology rebuild entry (not proof topology changed)',
 13:'scalar capture values produced: generation, delta physics values, delta geometry nodes, reserved',
 15:'RK tendency interval starts (stage marker precedes scalar production)',
 14:'catalog outcome: generation, total blocks, source-to-final owner changes, reserved',
 40:'selected collective site: site ID, datatype elements sent, received, reserved',
 41:'payload checksummed: generation, send bytes scanned, receive bytes scanned, reserved'}

def args_of(s):
    body=s[s.index('(')+1:s.rfind(')')];out=[];start=0;depth=0
    for i,c in enumerate(body):
        if c=='(':depth+=1
        if c==')':depth-=1
        if c==',' and depth==0:out.append(body[start:i].strip());start=i+1
    return out+[body[start:].strip()]

def event(kind,obj=0,gen='block_writeback_plan_generation',a='0',b='0',c='0'):
    return f'    if (lifecycle_enabled205()) call lifecycle_event205({kind},{obj}, &\n         int({gen},int64), &\n         int({a},int64), &\n         int({b},int64),int({c},int64))\n'

def instrument(root):
    root=Path(root);inventory={};site=0
    for p in sorted((root/'src').glob('*.f90')):
        if p.is_symlink() or p.name.startswith(('workflow_','parallel_block_profile')):continue
        original=p.read_text();lines=original.splitlines(True);out=[];edits=[];sites=[];i=0
        while i<len(lines):
            start=i;i+=1
            while lines[i-1].split('!')[0].rstrip().endswith('&') and i<len(lines):i+=1
            old=''.join(lines[start:i]);logical=' '.join(l.split('!')[0].strip().strip('&') for l in lines[start:i])
            m=re.match(r'\s*(?:if\s*\((.*)\)\s*)?call\s+(\w+)\b',logical,re.I)
            if not m or m[2].lower() not in TARGETS:out.append(old);continue
            target=m[2].lower();tag=TARGETS[target];indent=re.match(r'[ \t]*',old)[0];conditional=m[1]
            stmt=old
            if conditional is not None:
                # Preserve the entire original statement in restoration; evaluate
                # the original condition exactly once, before entering timers.
                pos=re.search(r'\bcall\s+',old,re.I).start();stmt=indent+old[pos:]
            begin=indent+f'call measure_enter205({tag})\n';end=indent+f'call measure_leave205({tag})\n'
            detail=''
            if tag in (26,27):
                site+=1;send=recv='-1';parts=args_of(logical[logical.lower().index('call'):]) if '(' in logical[logical.lower().index('call'):] else []
                if target=='mpi_alltoallv':send=f'sum({parts[1]})';recv=f'sum({parts[5]})'
                elif target=='mpi_allreduce':send=recv=parts[2]
                # Units are the site's MPI datatype elements; -1 means not counted.
                detail=f'{indent}if (lifecycle_selected205()) call lifecycle_event205(40,{site}, &\n{indent}     0_int64,int({send},int64), &\n{indent}     int({recv},int64),0_int64)\n'
            new=begin+detail+stmt+end
            if conditional is not None:new=indent+'if ('+conditional+') then\n'+new+indent+'end if\n'
            edits.append((new,old));out.append(new);sites.append({'line':start+1,'target':target,'operation':tag,'collective_site':site if tag in (26,27) else None,'statement':logical})
        changed=''.join(out)
        def replace_once(old,new):
            nonlocal changed
            if changed.count(old)!=1:raise ValueError(f'Anchor not unique {p.name}: {old[:100]!r}')
            changed=changed.replace(old,new);edits.append((new,old))
        if p.name=='parallel_block_mpi.f90':
            # Observe existing validity predicates. Equality of a plan generation
            # does not imply equality of prognostic fields or unchanged topology.
            for obj,old,cached in [(1,'      if (block_scalar_divergence_plan%generation == &\n           block_writeback_plan_generation) return','block_scalar_divergence_plan%generation'),
                 (2,'    if (block_scalar_restriction_exchange%generation == &\n         block_writeback_plan_generation) return','block_scalar_restriction_exchange%generation'),
                 (3,'    if (temperature_edge_generation == block_writeback_plan_generation) return','temperature_edge_generation')]:
                replace_once(old,event(1,obj,a=cached,b=f'merge(1,0,{cached}==block_writeback_plan_generation)')+old)
            old='    if (block_writeback_plan%ready) then\n       if (block_writeback_plan%catalog_size'
            replace_once(old,event(1,4,a='block_writeback_plan_generation',b='merge(1,0,block_writeback_plan%ready)')+old)
            old='    block_writeback_plan_generation = &\n         block_writeback_plan_generation + 1_int64'
            replace_once(old,old+'\n'+event(2,4))
            # clear is called both to invalidate and before initial construction;
            # ready distinguishes destruction of a valid plan from an empty clear.
            start=changed.index('  subroutine clear_block_writeback_plan');stop=changed.index('  end subroutine clear_block_writeback_plan',start)
            body=changed[start:stop];anchor=re.search(r'^    if \(allocated',body,re.M).start();old=body[anchor:].splitlines(True)[0]
            # Limit replacement to the routine because deallocation text recurs.
            newbody=body[:anchor]+event(3,4,a='merge(1,0,block_writeback_plan%ready)')+body[anchor:]
            replace_once(body,newbody)
            old='    if (block_scalar_divergence_plan%full_transport) then' # many sites: use exact setup anchor below
            anchor='    block_scalar_divergence_plan%producer_buffer = 0.0_dp\n'
            replace_once(anchor,anchor+event(4,1,a='merge(1,0,block_scalar_divergence_plan%full_transport)',b='merge(1,0,block_scalar_divergence_plan%production_cache_ready)',c='merge(1,0,storage_rebuilt)'))
            # Existing cumulative producer counters avoid adding point-level work.
            anchor='    call detail_enter(DP_SCALAR_SETUP)\n'
            replace_once(anchor,anchor+'    if (lifecycle_enabled205()) call lifecycle_producer_begin205(scalar_producer_work(2),scalar_producer_work(1))\n')
            anchor='    call detail_enter(DP_SCALAR_TRANSPORT)\n'
            replace_once(anchor,anchor+'    if (lifecycle_enabled205()) call lifecycle_producer_end205( &\n         block_writeback_plan_generation,scalar_producer_work(2),scalar_producer_work(1))\n')
            old='         block_scalar_divergence_plan%production_cache_ready = .true.'
            replace_once(old,old+'\n'+event(10,1))
            # Routine-local batch events, before numerical executable anchors.
            snippets=[('exchange_domain_to_block_payloads','    if (payload_family /= BLOCK_PAYLOAD_SOL',event(6,a='payload_family',b='merge(1,0,transfer_vector_only)')),
              ('exchange_block_writeback_payloads','    if (component == BLOCK_WRITEBACK_TEMPERATURE)',event(7,a='family',b='component')),
              ('capture_block_domain_multistage_candidate_tendency','    call block_profile_enter(BLOCK_PROFILE_TENDENCY)',event(8,a='stage',b='stage_count')),
              ('begin_block_domain_multistage_candidate_stage','    call detail_enter(DP_RK_ASSEMBLE)',event(9,a='stage',b='stage_count')),
              ('seed_scalar_producer_geometry','    call get_block_field_layout',event(5,1)),
              ('evaluate_candidate_block_scalar_restriction','    call detail_enter',event(11,1,a='level_start',b='level_end')),
              ('build_parallel_block_catalog','  call detail_enter',event(12))]
            for name,anchor,obs in snippets:
                start=changed.lower().index('  subroutine '+name);stop=changed.lower().index('end subroutine '+name,start)
                body=changed[start:stop];pos=body.index(anchor);line=body[pos:].splitlines(True)[0]
                replace_once(body,body[:pos]+obs+body[pos:])
            old='  n_changed = count(block_owner_new /= block_catalog%owner)'
            replace_once(old,old+'\n'+event(14,a='n_block_total',b='n_changed'))
            # Only the two full payload scans, excluding following reductions.
            start=changed.index('  subroutine exchange_block_migration_payloads');stop=changed.index('  end subroutine exchange_block_migration_payloads',start)
            body=changed[start:stop];a=body.index('    local_send_checksum = 0_int64');b=body.index('    call measure_enter205(27)',a)
            newbody=body[:a]+'    call measure_enter205(7)\n'+event(41,a='n_send_byte',b='n_recv_byte')+body[a:b]+'    call measure_leave205(7)\n'+body[b:]
            replace_once(body,newbody)
        if p.name=='time_integr.f90' and 'subroutine prepare_block_multistage_tendency' in changed:
            start=changed.index('  subroutine prepare_block_multistage_tendency')
            stop=changed.index('end subroutine prepare_block_multistage_tendency',start)
            body=changed[start:stop];pos=body.index('    call wf_enter(2)')
            replace_once(body,body[:pos]+event(15,gen='0',a='stage',b='stage_count')+body[pos:])
        if not sites and changed==original:continue
        m=re.search(r'^module \w+[^\n]*\n',changed,re.M|re.I)
        if not m:raise ValueError('Missing module import '+p.name)
        imp=('  use workflow_profile_195_mod, only: measure_enter205,measure_leave205,lifecycle_event205, &\n'
             '       lifecycle_enabled205,lifecycle_selected205,lifecycle_producer_begin205,lifecycle_producer_end205\n')
        # int64 is already imported in block files, but legacy MPI units need it.
        if not re.search(r'\b(?:use\s+iso_fortran_env[^\n]*int64|integer\s*,\s*parameter\s*::\s*int64)',original,re.I):
            imp+='  use iso_fortran_env, only: int64\n'
        new=m.group()+imp;edits.append((new,m.group()));changed=changed[:m.start()]+new+changed[m.end():]
        restored=changed
        for new,old in reversed(edits):restored=restored.replace(new,old)
        if restored!=original:raise ValueError('Numerical restoration failed '+p.name)
        p.write_text(changed);inventory[str(p.relative_to(root))]={'sites':sites,'original_bytes_recoverable':True}
    return inventory
