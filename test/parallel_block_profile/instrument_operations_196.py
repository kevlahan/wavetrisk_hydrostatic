"""Append observer-only calls; every original source byte must survive removal."""
from pathlib import Path
import re
from prepare_legacy_profile import replace_once

# These wrappers gate old detail hooks differently. Emit independent, paired
# hooks before their old profiling guard; leave existing code unaltered.
WRAPPERS=('parallel_block_profile_begin','parallel_block_profile_end','block_profile_enter','block_profile_leave')

def instrument(root,block):
    root=Path(root);changed=[];inventory={}
    for p in sorted((root/'src').glob('*.f90')):
        if p.is_symlink() or p.name=='workflow_profile_195.f90':continue
        before=p.read_text();s=before;additions=[];counts={}
        def add(old,new):
            nonlocal s
            s=replace_once(s,old,new);additions.append((new,old))
        if block and p.name=='parallel_block_mpi.f90':
            protected=[]
            for name in WRAPPERS:
                match=re.search(r'^  (?:real\(dp\) function|subroutine) '+name+r'\b.*?^  end (?:function|subroutine) '+name+r'\b',s,re.M|re.S)
                if not match:raise ValueError('Missing wrapper '+name)
                body=match.group();case=re.search(r'    select case\(phase\).*?    end select',body,re.S).group()
                # Exclude the reset guard in the step case; only route observer calls.
                proxy='    select case(phase)\n'
                chunks=re.split(r'(?m)^    case',case)[1:]
                for chunk in chunks:
                    selector=chunk.splitlines()[0]
                    hit=re.search(r'call detail_(enter|leave)\(([^\n]+)\)',chunk)
                    if not hit:raise ValueError('Unrecognized wrapper case')
                    proxy+='    case'+selector+'\n       call op_'+hit[1]+'('+hit[2]+')\n'
                proxy+='    end select\n'
                anchor='    if (phase < 1 .or. phase > BLOCK_PROFILE_PHASE_COUNT)'
                new=replace_once(body,anchor,proxy+anchor)
                add(body,new)
                protected.append((name,new))
            # Direct detail scopes outside the wrappers get their own operation hook.
            placeholders={}
            for i,(name,body) in enumerate(protected):
                key=f'! OPERATIONS196_PROTECTED_{i}\n';s=s.replace(body,key);placeholders[key]=body
        else:placeholders={}
        # Calls are whole single-line statements in the frozen sources.
        pattern=r'^(\s*)call detail_(enter|leave)\(([^\n]+)\)[ \t]*$'
        def hook(m):
            name=m[2];region=m[3];counts[region]=counts.get(region,0)+1
            # Preserve original indentation/newlines verbatim in the original statement.
            indent=m[1].split('\n')[-1]
            return m.group()+f'\n{indent}call op_{name}({region})'
        s=re.sub(pattern,hook,s,flags=re.M)
        for key,body in placeholders.items():s=s.replace(key,body)
        if p.name=='parallel_block_profile.f90':
            # Work counters are already computed at batch/route granularity.
            a='    if (p<1.or.p>DC_COUNT.or.n<0)';idx=s.index('  subroutine detail_add');tail=s[idx:]
            target=tail[:tail.index('  end subroutine')]
            new=replace_once(target,'    if (.not.detail_enabled) return','    call op_count(p,n)\n    if (.not.detail_enabled) return')
            add(target,new)
        if block and p.name=='parallel_block_mpi.f90':
            idx=s.index('  subroutine record_parallel_block_profile_volume');end=s.index('  end subroutine record_parallel_block_profile_volume',idx)
            body=s[idx:end];new=replace_once(body,'    if (.not. parallel_block_profile_enabled()) return',
                 '    call op_volume(phase,message_count,byte_count)\n    if (.not. parallel_block_profile_enabled()) return')
            add(body,new)
        if block and p.name=='parallel_block_mpi.f90':
            # The old detail profiler did not separate these reverse imports.
            # Instrument complete outer bodies; internal recursive routines
            # remain untouched and no per-patch clock reads are introduced.
            for name,first,id in [
                ('import_domain_field_family_to_blocks','    keep_checkpoint = .false.',70),
                ('exchange_domain_to_block_payloads','    use_tags = .false.',71),
                ('import_domain_boundary_field_family_to_blocks','    check_current = .false.',73),
                ('write_block_field_family_to_domains','    call block_profile_enter(BLOCK_PROFILE_WRITEBACK)',77)]:
                idx=s.index('  subroutine '+name+' (');end=s.index('  end subroutine '+name,idx)
                body=s[idx:end]
                stop=body.index('  contains') if '  contains' in body else len(body)
                if re.search(r'\breturn\b', '\n'.join(l for l in body[:stop].splitlines() if not l.lstrip().startswith('!'))):
                    raise ValueError('Unexpected early return in '+name)
                head=replace_once(body[:stop],first,f'    call op_enter({id})\n'+first)
                if id==71:
                    anchor='    if (.not. transfer_vector_only) &\n         block_writeback_plan%scalar_recv_buffer = 0.0_dp'
                    head=replace_once(head,anchor,'''    call op_count(23, &
         merge(int(size(block_writeback_plan%scalar_recv_buffer),int64),0_int64,.not.transfer_vector_only) + &
         merge(int(size(block_writeback_plan%scalar_send_buffer),int64),0_int64,.not.transfer_vector_only) + &
         int(size(block_writeback_plan%vector_recv_buffer),int64) + &
         int(size(block_writeback_plan%vector_send_buffer),int64))
'''+anchor)
                if id==70:
                    head += '    call op_count(24,int(imported_patch_count,int64))\n'
                if id==73:
                    anchor='    if (transfer_scalar) then\n       block_writeback_plan%boundary_scalar_domain_send_buffer = 0.0_dp'
                    head=replace_once(head,anchor,'''    call op_leave(73)
    call op_enter(74)
    call op_count(25, &
         merge(int(size(block_writeback_plan%boundary_scalar_domain_send_buffer),int64) + &
               int(size(block_writeback_plan%boundary_scalar_block_recv_buffer),int64),0_int64,transfer_scalar) + &
         merge(int(size(block_writeback_plan%boundary_vector_domain_send_buffer),int64) + &
               int(size(block_writeback_plan%boundary_vector_block_recv_buffer),int64),0_int64,transfer_vector))
'''+anchor)
                    anchor='    do r = 1,n_process\n       pos_scalar = block_writeback_plan% &\n            boundary_scalar_block_recv_displ(r) + 1'
                    head=replace_once(head,anchor,'    call op_leave(74)\n    call op_enter(76)\n'+anchor)
                    anchor='       n_boundary = local_block_boundary_count(b)\n       do boundary_index = 1,n_boundary'
                    head=replace_once(head,anchor,'       n_boundary = local_block_boundary_count(b)\n       call op_count(26,int(n_boundary,int64))\n       do boundary_index = 1,n_boundary')
                    anchor='          if (n_boundary /= &\n               block_writeback_plan%send_boundary_count(slot)) then'
                    head=replace_once(head,anchor,'          call op_count(27,int(n_boundary,int64))\n'+anchor)
                if id in (71,73):
                    mpi_id=72 if id==71 else 75
                    pat=r'^([ \t]*)call MPI_Alltoallv\( &\n.*?MPI_DOUBLE_PRECISION,comm,ierr\)'
                    def mpi_hook(m):
                        return m[1]+f'call op_enter({mpi_id})\n'+m.group()+'\n'+m[1]+f'call op_leave({mpi_id})'
                    head,nmpi=re.subn(pat,mpi_hook,head,flags=re.M|re.S)
                    if nmpi!=2:raise ValueError('Expected two payload collectives in '+name)
                head+=f'    call op_leave({76 if id==73 else id})\n'
                add(body,head+body[stop:])
        if block and p.name=='multi_level.f90':
            for name,id in [('begin_native_mass',68),('finish_native_mass',69)]:
                idx=s.index('  subroutine '+name+'(');end=s.index('  end subroutine',idx)
                body=s[idx:end];new=replace_once(body,'    associate(p=>mass_boundary(pos,l))',f'    call op_enter({id})\n    associate(p=>mass_boundary(pos,l))')
                new+=f'    call op_leave({id})\n'
                if id==68:
                    anchor='    p%nrequest=nreq'
                    new=replace_once(new,anchor,'''    call op_count(19,int(size(p%send_value),int64))
    call op_count(20,int(size(p%recv_value),int64))
    call op_count(21,int(size(p%local,2),int64))
    call op_count(22,int(nreq,int64))
'''+anchor)
                add(body,new)
        if s!=before:
            if 'use workflow_profile_195_mod' not in s:
                m=re.search(r'^module \w+[^\n]*\n',s,re.M)
                if not m:raise ValueError('Cannot import observer in '+p.name)
                anchor=s[m.start():m.end()+80]
                add(anchor,m.group()+'  use workflow_profile_195_mod, only: op_enter,op_leave,op_count,op_volume\n'+s[m.end():m.end()+80])
            # Remove exact structural additions first, then simple appended hooks.
            restored=s
            for new,old in reversed(additions):restored=replace_once(restored,new,old)
            restored=''.join(l for l in restored.splitlines(True) if not re.match(r'^\s*call op_(enter|leave)\(',l))
            if restored!=before:raise ValueError('Numerical source changed: '+p.name)
            p.write_text(s);changed.append('src/'+p.name);inventory[p.name]=counts
    (root/'src/workflow_profile_195.f90').write_bytes(Path(__file__).with_name('workflow_profile_196.f90').read_bytes())
    return {'changed_sources':changed,'direct_scope_hooks':inventory,'scope':'Observer additions only; stripping restores every original numerical source byte.'}
