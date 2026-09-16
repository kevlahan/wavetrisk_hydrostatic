"""Batch-level causal diagnostics over frozen Stage198 observer sources."""
from pathlib import Path
import re
from prepare_legacy_profile import replace_once
REGIONS={
101:'basic operators batch',102:'pressure and step1 patch batch',103:'scalar physics capture patch batch',
104:'scalar divergence capture batch',105:'Bernoulli restriction batch',106:'Exner restriction batch',
107:'native mass restriction batch',108:'mass exchange post and pack',109:'mass exchange wait and install',
110:'scalar direct flux all-block batch',111:'scalar level restriction batch',112:'scalar level divergence batch',
113:'temperature edge exchange batch',114:'scalar ghost exchange batch',115:'scalar capture transport batch',
116:'catalog local traversal',117:'catalog descriptor allocation',118:'source store allocation',
119:'source extraction batch',120:'migration route construction',121:'migration sizing traversal',
122:'migration buffer allocation',123:'migration serialization batch',124:'migration receive allocation',
125:'migration local installation',126:'migration retained-state preparation',127:'migration staging release',
128:'source field/topology allocation',129:'source topology and route construction',130:'scalar mass RHS publication',131:'source interior field copying and validation',
132:'source boundary construction and field copying',133:'source ghost construction and field copying',
134:'source stencil construction and validation',135:'catalog ownership balancing'}
COUNTERS={}

def upgrade(root,block):
 root=Path(root);inventory={}
 for p in sorted((root/'src').glob('*.f90')):
  if p.is_symlink() or p.name=='workflow_profile_195.f90':continue
  before=p.read_text();s=before;edits=[];counts={}
  def change(old,new):
   nonlocal s
   if old not in s:raise ValueError('Missing Stage199 anchor in '+p.name+': '+old[:90])
   if not new:new='! stage199 removed observer '+str(len(edits))+'\n'
   s=s.replace(old,new);edits.append((new,old))
  def wrap(old,id,args=''):
   indent=re.match(r'[ \t]*',old)[0]
   new=indent+f'call batch_begin199({id}{args})\n'+old+('\n' if not old.endswith('\n') else '')+indent+f'call batch_end199({id})\n'
   change(old,new);counts[id]=counts.get(id,0)+1
  def routine(name):
   m=re.search(r'^  subroutine '+name+r'\b.*?^  end subroutine '+name+r'\b[^\n]*',s,re.M|re.S)
   if not m:raise ValueError('Missing routine '+name)
   return m.group()
  def calls(text,pattern,id,args=''):
   # Match entire call, including continuations and optional one-line condition.
   lines=text.splitlines(True);i=0;seen=set()
   while i<len(lines):
    start=i;i+=1
    if not re.search(pattern,lines[start],re.I):continue
    while lines[i-1].split('!')[0].rstrip().endswith('&'):i+=1
    old=''.join(lines[start:i])
    if old not in seen:wrap(old,id,args);seen.add(old)
  if p.name=='multi_level.f90':
   old=routine('basic_operators');body=old
   body=replace_once(body,'    do d = 1, size(grid)','    call batch_begin199(101,k,l)\n    do d = 1, size(grid)')
   body=body.replace('  end subroutine basic_operators','    call batch_end199(101)\n  end subroutine basic_operators')
   if block:
    body=body.replace('    call batch_end199(101)\n  end subroutine basic_operators','  end subroutine basic_operators')
    body=replace_once(body,'    call detail_leave(DP_BASIC)','    call batch_end199(101)\n    call detail_leave(DP_BASIC)')
   change(old,body)
   # Existing per-patch legacy step1 observers would otherwise hide the matched
   # producer under a different label. The identical call bytes are preserved.
   body=routine('basic_operators')
   for line in ('          call op_enter(79)\n','          call op_leave(79)\n'):
    if line in body:change(body,body.replace(line,''));body=routine('basic_operators')
   calls(body,r'call apply_onescale_to_patch \(integrate_pressure_up',102,',k,l,work=PATCH_SIZE**2')
   calls(body,r'call step1\s*\(',102,',k,l,work=PATCH_SIZE**2')
   calls(body,r'call capture_block_scalar_physics_patch',103,',k,l,work=PATCH_SIZE**2')
   # Scalar restriction is two separate producers even though it calls one abstraction.
   body=routine('basic_operators')
   for field,id in [('grid(d)%bernoulli%elts',105),('exner_fun(k)%data(d)%elts',106)]:
    a=body.index('          scalar => '+field);b=body.index('          nullify (scalar)',a)+len('          nullify (scalar)')
    wrap(body[a:b],id,',k,l,work=grid(d)%lev(l)%length');body=routine('basic_operators')
   if block:
    calls(s,r'call capture_block_scalar_divergence_level',104,',k,l')
    body=routine('basic_operators');a=body.index('       do d=1,size(native_mass)');b=body.index('       if (level_start/=level_end)',a)
    wrap(body[a:b],107,',k,l,work=size(native_mass)')
    # The begin/finish routines have no early return. Their outer batches carry
    # exact source/local/received counts and preserve async posting order.
    old=routine('begin_native_mass');new=old.replace('    if(p%pending)', '    call batch_begin199(108,k,l,pos,size(p%send_value)+size(p%local,2))\n    if(p%pending)',1)
    new=new.replace('    end associate','    call batch_end199(108)\n    end associate',1);change(old,new)
    # finish_native_mass has an unnamed END; locate by its next routine.
    a=s.index('  subroutine finish_native_mass(');b=s.index('  subroutine exchange_native_mass(',a)
    old=s[a:b];new=old.replace('    if(.not.p%pending','    call batch_begin199(109,k,l,pos,size(p%recv_value))\n    if(.not.p%pending',1)
    new=new.replace('    end associate','    call batch_end199(109)\n    end associate',1);change(old,new)
  if block and p.name=='parallel_block_mpi.f90':
   # All-layer batches use k=0, unknown/all-level keys use l=0 or -1.
   body=routine('evaluate_candidate_block_scalar_restriction')
   for target,id,args in [('compute_block_scalar_direct_flux',110,',0,level_end'),('restrict_block_scalar_flux_level',111,',0,statistics%target_level'),('recompute_block_scalar_divergence_level',112,',0,statistics%target_level')]:
    pat=r'(?m)^([ \t]*)call apply_local_block_field_consumer\( &\n[ \t]*'+target+r'[^\n]*\n'
    for m in list(re.finditer(pat,body)):wrap(m.group(),id,args)
   for lev in ('level_end','l'):
    for m in list(re.finditer(r'(?m)^([ \t]*)call exchange_temperature_edges\('+lev+r'\)[^\n]*\n',s)):
     wrap(m.group(),113,',0,'+lev)
   body=routine('exchange_block_scalar_restriction_ghosts')
   anchor='    call block_profile_enter(profile_ghost_mode)'
   new=replace_once(body,anchor,'    call batch_begin199(114,0,exchange_level,component)\n'+anchor)
   new=replace_once(new,'    call detail_leave(DP_GHOST)','    call batch_end199(114)\n    call detail_leave(DP_GHOST)');change(body,new)
   calls(s,r'call publish_native_mass_rhs\b',130,',0,0')
   body=s[s.index('  subroutine build_parallel_block_catalog'):s.index('end subroutine build_parallel_block_catalog')]
   a=body.index('  n_block_local = 0');b=body.index('  allocate(block_loc(n_block_local))',a)
   wrap(body[a:b],116);wrap('  allocate(block_loc(n_block_local))\n',117,',work=n_block_local')
   body=s[s.index('  subroutine build_parallel_block_catalog'):s.index('end subroutine build_parallel_block_catalog')]
   a=body.index('  i = 0');b=body.index('  call MPI_Exscan',a)
   # Exclude any preceding existing MPI scope opener from the traversal.
   tail=body[a:b];cut=tail.rfind('  call op_enter(')
   if cut>=0:tail=tail[:cut]
   wrap(tail,116)
   a=s.index('  block_owner_new = -1');b=s.index('  detailed_diagnostics = block_detailed_diagnostics_enabled()',a)
   wrap(s[a:b],135)
   body=routine('migrate_blocks')
   for name,id in [('build_block_migration_manifest',120),('install_local_blocks',125),('prepare_parallel_block_production_state',126),('clear_block_staging',127),('clear_block_migration_manifest',127)]:calls(body,r'call '+name+r'\b',id)
   body=routine('migrate_blocks');a=body.index('    allocate(send_nbyte(');b=body.index('    call exchange_block_migration_sizes',a)
   # Existing op scopes before exchange belong to MPI, not sizing.
   chunk=body[a:b];cut=chunk.rfind('    call op_enter(')
   if cut>=0:chunk=chunk[:cut]
   wrap(chunk,121)
   body=routine('migrate_blocks');a=body.index('    allocate(send_payload(');b=body.index('    pos = 0',a);wrap(body[a:b],122)
   body=routine('migrate_blocks');a=body.index('    do i = 1, manifest%n_send',body.index('    allocate(send_payload('));b=body.index('    call detail_leave(DP_PACK)',a);wrap(body[a:b],123)
   body=routine('migrate_blocks');a=body.index('    if (allocated(block_received))');b=body.index('    pos = 0',a);wrap(body[a:b],124)
  if block and p.name=='parallel_block_build.f90':
   a=s.index('  if (allocated(block_source))');b=s.index('  b_verbose = -1',a);wrap(s[a:b],118)
   calls(s,r'call build_one_source_block\b',119,',work=block_catalog(b)%weight')
   a=s.index('  call get_block_field_layout(');b=s.index('  ! Verify child links, leaf count and depth.',a)
   wrap(s[a:b],131)
   a=s.index('  allocate(block_out%neigh_class(');b=s.index('  allocate(ghost_patch(',a)
   wrap(s[a:b],132)
   a=s.index('  allocate(ghost_patch(');b=s.index('  call renumber_subtree_neigh_Domain(',a)
   wrap(s[a:b],133)
   a=s.index('  allocate(block_out%stencil(');b=s.index('\ncontains\n',a)
   wrap(s[a:b],134)
   # Contiguous allocation statements, including continuation lines. Never
   # wrap allocation predicates or declaration attributes.
   calls(s,r'^\s*allocate\(',128)
   for name in ('extract_subtree_patches_Domain','compact_subtree_storage_Domain','renumber_subtree_neigh_Domain'):
    calls(s,r'call '+name+r'\b',129)
  if block:
   calls(s,r'call finalize_block_scalar_divergence_capture\b',115,',0,0')
  if block and p.name=='parallel_block_scalar_storage.f90':
   # Keep only guarded integer census. No clocks or sampler calls inside reads/writes.
   for m in list(re.finditer(r'#ifdef WAVETRISK_STORAGE198\n.*?#endif\n',s,re.S)):
    old=m.group()
    if old not in s:continue
    kept=''.join(line for line in old.splitlines(True) if 'if (count198)' in line or 'if (count198.and.' in line)
    change(old,('#ifdef WAVETRISK_STORAGE199\n'+kept.replace('198','199').replace('count199','census199')+'#endif\n') if kept else '')
  if s!=before:
   if 'use workflow_profile_195_mod' not in s:
    m=re.search(r'^module \w+[^\n]*\n',s,re.M);assert m;change(m.group(),m.group()+'  use workflow_profile_195_mod\n')
   for m in list(re.finditer(r'(?m)^  use workflow_profile_195_mod, only:[^\n]*',s)):change(m.group(),'  use workflow_profile_195_mod')
  # The preceding observer names change only; numerical source bytes are unchanged.
  if 'count198' in s or 'leaf198' in s:change(s,s.replace('count198','count199').replace('leaf198','leaf199'))
  restored=s
  for new,old in reversed(edits):restored=restored.replace(new,old)
  if restored!=before:raise ValueError('Stage199 hooks cannot restore Stage198 bytes: '+p.name)
  if s!=before:p.write_text(s);inventory[p.name]={'batch_sites':counts,'restores_stage198_bytes':True}
 make=root/'Makefile';make.write_text(make.read_text().replace('workflow_clock_198','workflow_clock_199'))
 (root/'src/workflow_clock_198.c').unlink()
 (root/'src/workflow_clock_199.c').write_bytes(Path(__file__).with_name('workflow_clock_199.c').read_bytes())
 (root/'src/workflow_profile_195.f90').write_bytes(Path(__file__).with_name('workflow_profile_199.f90').read_bytes())
 return inventory
