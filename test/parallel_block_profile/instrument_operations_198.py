"""Diagnosis-only additions to frozen Stage197 observers; no numerical edits."""
from pathlib import Path
import re
from prepare_legacy_profile import replace_once
REGIONS={83:'catalog construction',84:'source-block construction',85:'block migration/install',86:'post-grid reconstruction retention',87:'adaptation lifecycle completion',88:'Domain inverse transform envelope',89:'mask/threshold updates',90:'surface pressure producer',91:'native velocity plan preparation',92:'native mass plan preparation',93:'scalar address/geometry preparation',94:'scalar geometry seeding',95:'scalar capture finalization',96:'complete vector wavelet producer',97:'compression entry envelope',98:'native RK retention',99:'compatibility RK assembly',100:'diffusion preparatory operators'}
COUNTERS={36:'scalar range requests',37:'scalar range logical slots requested',38:'complete 50-slot ranges requested',39:'scalar read_value visits',40:'shared geometry address resolutions',41:'scalar write_value visits',42:'scalar column read requests'}
TARGETS={83:['build_parallel_block_catalog'],84:['build_source_blocks'],85:['migrate_blocks'],86:['retain_post_grid_change_block_reconstruction'],87:['complete_block_adaptation_lifecycle'],88:['inverse_wavelet_transform'],89:['init_masks_zero','mask_active','mask_adj_same_scale','mask_restrict_same_scale','mask_adj_finer_scale','complete_masks','mask_second_neighbours'],90:['cal_surf_press'],91:['prepare_native_velocity_programs'],92:['prepare_native_mass'],93:['prepare_block_scalar_capture_addresses'],94:['seed_scalar_producer_geometry'],95:['finalize_block_scalar_divergence_capture'],96:['validate_complete_block_vector_wavelets'],97:['compress_wavelets'],98:['retain_block_native_multistage_candidate'],99:['RK_sub_step_compatibility'],100:['cal_divu_ml','cal_Laplacian_scalars','cal_Laplacian_divu']}

def upgrade(root,block):
 root=Path(root);inventory={};wanted={name.lower():i for i,names in TARGETS.items() for name in names}
 for p in sorted((root/'src').glob('*.f90')):
  if p.is_symlink() or p.name=='workflow_profile_195.f90':continue
  before=p.read_text();s=before;changes=[];counts={}
  def change(old,new):
   nonlocal s
   pattern=re.compile('^'+re.escape(old),re.M)
   if not pattern.search(s):raise ValueError('Missing observer anchor '+old[:100])
   s=pattern.sub(lambda m:new,s);changes.append((new,old))
  # Match complete Fortran call statements, preserving original bytes and conditionals.
  lines=s.splitlines(True);i=0;seen_calls=set()
  while i<len(lines):
   first=i;line=lines[i];m=re.match(r'^([ \t]*)(?:if\s*\([^\n]*\)\s*)?call\s+(\w+)\b',line,re.I)
   i+=1
   if not m or m[2].lower() not in wanted:continue
   while lines[i-1].split('!')[0].rstrip().endswith('&') and i<len(lines):i+=1
   old=''.join(lines[first:i]);indent=m[1];region=wanted[m[2].lower()]
   if old in seen_calls:continue
   seen_calls.add(old)
   counts[region]=counts.get(region,0)+len(re.findall('^'+re.escape(old),s,re.M))
   new=indent+f'call op_enter({region})\n'+old+indent+f'call op_leave({region})\n'
   change(old,new)
  if not block and p.name=='ops.f90':
   # Plain integer accumulation replaces per-leaf calls; flush at common context boundaries.
   pattern=r'(?m)^([ \t]*)call op_count\((3[1345]),([^\n]+)\)$'
   seen_leaf=set()
   for m in list(re.finditer(pattern,s)):
    old=m.group()
    if old in seen_leaf:continue
    seen_leaf.add(old)
    new=m[1]+f'if (count198) leaf198({m[2]})=leaf198({m[2]})+{m[3]}'
    change(old,new)
  if block and p.name=='parallel_block_scalar_storage.f90':
   for name,probe,counter,anchor in [
    ('scalar_read_range',1,36,'    if (first<1.or.last>SCALAR_RECORD_WIDTH*store%samples)'),
    ('read_value',2,39,'    if (address<1.or.address>SCALAR_RECORD_WIDTH*store%samples)'),
    ('shared_address',3,40,'    field=mod(sample/store%stride,store%nfield)'),
    ('write_value',4,41,'    if (address<1.or.address>SCALAR_RECORD_WIDTH*store%samples)')]:
    match=re.search(r'^  (?:[^\n]*function|subroutine) '+name+r'\(.*?^  end (?:function|subroutine)(?: '+name+r')?[^\n]*',s,re.M|re.S)
    if not match:raise ValueError('Missing storage routine '+name)
    old=match.group();body=old
    # Every return must close a sampled invocation; no clocks in normal observer builds.
    start=body.index(anchor) if anchor in body else -1
    if start<0:raise ValueError('Missing storage anchor '+name)
    intro='''#ifdef WAVETRISK_STORAGE198
    logical :: take198
    take198=sample_begin198(PROBE)
    if (count198) leaf198(COUNT)=leaf198(COUNT)+1_int64
'''.replace('PROBE',str(probe)).replace('COUNT',str(counter))
    if name=='scalar_read_range':intro+='    if (count198) leaf198(37)=leaf198(37)+int(max(0,last-first+1),int64)\n    if (count198.and.last-first+1==50) leaf198(38)=leaf198(38)+1_int64\n'
    intro+='#endif\n'
    body=body[:start]+intro+body[start:]
    ret=r'(?m)^([ \t]*)(?:if\s*(\([^\n]*\))\s*)?return\s*$'
    def close(m):
     indent,cond=m[1],m[2]
     if cond:return indent+'if '+cond+' then\n#ifdef WAVETRISK_STORAGE198\n'+indent+f'  call sample_end198({probe},take198)\n#endif\n'+indent+'  return\n'+indent+'end if'
     return '#ifdef WAVETRISK_STORAGE198\n'+indent+f'call sample_end198({probe},take198)\n#endif\n'+m.group()
    body=re.sub(ret,close,body)
    end=body.rindex('  end ');body=body[:end]+f'#ifdef WAVETRISK_STORAGE198\n    call sample_end198({probe},take198)\n#endif\n'+body[end:]
    change(old,body)
   a='    call column_bounds(store,sample,slot,count,last,mapped)'
   # A column read count only; write-column calls share this helper and are excluded.
   start=s.index('  function scalar_read_column(');end=s.index('  end function scalar_read_column',start)
   old=s[start:end];new=replace_once(old,a,'#ifdef WAVETRISK_STORAGE198\n    if (count198) leaf198(42)=leaf198(42)+1_int64\n#endif\n'+a);change(old,new)
  if s!=before:
   if 'use workflow_profile_195_mod' not in s:
    m=re.search(r'^module \w+[^\n]*\n',s,re.M);assert m
    change(m.group(),m.group()+'  use workflow_profile_195_mod\n')
   else:
    for m in list(re.finditer(r'(?m)^  use workflow_profile_195_mod, only:[^\n]*',s)):change(m.group(),'  use workflow_profile_195_mod')
   restored=s
   for new,old in reversed(changes):restored=re.sub('^'+re.escape(new),lambda m:old,restored,flags=re.M)
   if restored!=before:raise ValueError('Observer additions cannot restore original bytes: '+p.name)
   p.write_text(s);inventory[p.name]={'batch_scope_calls':counts,'restores_stage197_bytes':True}
 makefile=root/'Makefile';before=makefile.read_text()
 anchor='OBJ := $(addprefix $(BUILD_DIR)/,$(SRC:.f90=.o))'
 replacement=anchor+' $(BUILD_DIR)/workflow_clock_198.o'
 rule='\nMPICC ?= mpicc\n\n$(BUILD_DIR)/workflow_clock_198.o: $(SRC_DIR)/workflow_clock_198.c | dirs\n\t$(MPICC) -O2 -c $< -o $@\n'
 after=replace_once(before,anchor,replacement)+rule
 if after.removesuffix(rule).replace(replacement,anchor)!=before:raise ValueError('Clock build rules cannot restore original bytes')
 makefile.write_text(after)
 (root/'src/workflow_clock_198.c').write_bytes(Path(__file__).with_name('workflow_clock_198.c').read_bytes())
 inventory['Makefile']={'restores_stage197_bytes':True,'purpose':'identical monotonic observer clock'}
 (root/'src/workflow_profile_195.f90').write_bytes(Path(__file__).with_name('workflow_profile_198.f90').read_bytes())
 return inventory
