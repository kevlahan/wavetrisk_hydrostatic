"""Focused observers added to a Stage196-instrumented disposable source copy."""
from pathlib import Path
import re
from prepare_legacy_profile import replace_once
REGIONS={78:'resident coefficient initialization',79:'direct scalar flux and its enclosing traversal',80:'scalar or flux restriction traversal',81:'scalar divergence/tendency traversal',82:'temperature boundary closure'}
COUNTERS={28:'resident coefficient doubles zeroed',29:'scalar coefficient imports avoided',30:'uncompressed publications avoided',31:'direct scalar node-field samples',32:'restricted scalar flux edges',33:'scalar divergence samples',34:'scalar divergence evaluations',35:'temperature direct-flux node evaluations'}

def upgrade(root,block):
 root=Path(root);changes={}
 for p in sorted((root/'src').glob('*.f90')):
  if p.is_symlink() or p.name=='workflow_profile_195.f90':continue
  original=p.read_text();s=original
  # Call-site guards must independently enable the new counters.
  if p.name in ('comm_mpi.f90','parallel_block_mpi.f90'):
   s=s.replace('if (detail_enabled)', 'if (detail_enabled.or.op_enabled())')
   s=s.replace('if (detail_enabled.and.', 'if ((detail_enabled.or.op_enabled()).and.')
  def wrap_call(pattern,region):
   nonlocal s
   def wrap(m):
    indent=m[1];return indent+f'call op_enter({region})\n'+m.group()+'\n'+indent+f'call op_leave({region})'
   s,n=re.subn(pattern,wrap,s,flags=re.M)
   if not n:raise ValueError('Missing focused call '+pattern)
  if block and p.name=='parallel_block_mpi.f90':
   if 'subroutine zero_resident_wavelets' in s:
    wrap_call(r'^( +)call apply_local_block_field_producer\(zero_resident_wavelets,statistics\)',78)
    s=replace_once(s,'       block%wavelet_vector=0.0_dp', '       block%wavelet_vector=0.0_dp\n       call op_count(28,int(size(block%wavelet_scalar),int64)+int(size(block%wavelet_vector),int64))')
    s=replace_once(s,'    stage197_coefficients_ready=stage197_resident', '    if (stage197_resident.and..not.block_dynamics_validation_enabled()) call op_count(30,1_int64)\n    stage197_coefficients_ready=stage197_resident')
    s=replace_once(s,'    if (stage197_resident) then\n       if (.not.stage197_coefficients_ready', '    if (stage197_resident) then\n       call op_count(29,1_int64)\n       if (.not.stage197_coefficients_ready')
   for name,i in [('compute_block_scalar_direct_flux',79),('restrict_block_scalar_flux_level',80),('recompute_block_scalar_divergence_level',81),('compute_temperature_rk_boundary',82)]:
    wrap_call(r'^( +)call apply_local_block_field_consumer\((?: &\n +)?'+name+r',[^\n]+\)',i)
   a=s.index('  subroutine compute_block_scalar_direct_flux (');b=s.index('  end subroutine compute_block_scalar_direct_flux',a)
   t=s[a:b];t=replace_once(t,'                   i = mod(q,PATCH_SIZE)',
       '                   if (block%scalar_variable+scalar_slot==S_TEMP) &\n                        statistics%scalar_sample_count=statistics%scalar_sample_count+1_int64\n                   i = mod(q,PATCH_SIZE)')
   s=s[:a]+t+s[b:]
   a=s.index('  subroutine recompute_block_scalar_divergence_level (');b=s.index('  end subroutine recompute_block_scalar_divergence_level',a)
   t=s[a:b]
   t=replace_once(t,'    logical :: validate_oracle','    logical :: validate_oracle\n    integer(int64) :: measured_samples,measured_evaluations')
   t=replace_once(t,'    validate_oracle=block_dynamics_validation_enabled()', '    measured_samples=0_int64;measured_evaluations=0_int64\n    validate_oracle=block_dynamics_validation_enabled()')
   t=replace_once(t,'                   do k=1,count','                   measured_samples=measured_samples+int(count,int64)\n                   do k=1,count\n                      if (column(5,k)>0.5_dp) measured_evaluations=measured_evaluations+1_int64')
   t=replace_once(t,'    end select\n\n  contains','    end select\n    call op_count(33,measured_samples)\n    call op_count(34,measured_evaluations)\n\n  contains')
   s=s[:a]+t+s[b:]
   anchor='    count_local = [statistics%parent_patch_count, &\n         statistics%restricted_edge_count, &'
   s=replace_once(s,anchor,'''    call op_count(31,direct_flux_statistics%sample_count)
    call op_count(32,statistics%restricted_edge_count)
    call op_count(35,direct_flux_statistics%scalar_sample_count)
'''+anchor)
  if not block and p.name=='multi_level.f90':
   wrap_call(r'^( +)call cpt_or_restr_(?:scalar|flux) \([^\n]+\)[^\n]*',80)
   wrap_call(r'^( +)call apply_onescale_to_patch \(scalar_trend,[^\n]+\)',81)
   # step1 contains shared mass/velocity work too; do not equate its wall ratio to the native direct-flux loop.
   wrap_call(r'^( +)call step1 \([^\n]+\)',79)
   anchor='    id_par = idx (i_par, j_par, offs_par, dims_par)\n    \n    if (maxval(dom%mask_e'
   s=replace_once(s,anchor,'''    id_par = idx (i_par, j_par, offs_par, dims_par)
    call op_count(32,int(count(dom%mask_e%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)>=RESTRCT),int64))

    if (maxval(dom%mask_e''')
  if not block and p.name=='ops.f90':
   anchor='    if (dom%mask_n%elts(id_i) >= TRSK) then\n       dscalar(id_i) = - div'
   s=replace_once(s,anchor,'''    call op_count(33,1_int64)
    if (dom%mask_n%elts(id_i) >= TRSK) then
       call op_count(34,1_int64)
       dscalar(id_i) = - div''')
   anchor='         do v = scalars(1), scalars(2)\n            physics_flux = physics_scalar_flux'
   s=s.replace(anchor,'''         call op_count(35,1_int64)
         call op_count(31,int(scalars(2)-scalars(1)+1,int64))
         do v = scalars(1), scalars(2)
            physics_flux = physics_scalar_flux''')
  # One inverse exchange at a time. Dependency records include no-request completions.
  if block and p.name=='parallel_block_inverse.f90':
   s=replace_once(s,'    call detail_leave(DP_INVERSE_POST)','    call op_post(nrequest,n*sum(sc),n*sum(rc),tag)\n    call detail_leave(DP_INVERSE_POST)')
   s=replace_once(s,'    call detail_enter(DP_INVERSE_WAIT)','    call op_wait_begin\n    call detail_enter(DP_INVERSE_WAIT)')
   s=replace_once(s,'    call detail_leave(DP_INVERSE_WAIT)','    call detail_leave(DP_INVERSE_WAIT)\n    call op_wait_end')
   for name,first,kind in [('native_inverse_gather','    reset=.true.',1),('native_inverse_scatter','    if (present(changed_level)) then',2),('native_inverse_boundary','    call detail_enter',3)]:
    a=s.index('  subroutine '+name+'(');b=s.index('  end subroutine',a);t=s[a:b]
    if kind<3:
     insert=f'    call op_dependency({kind},component,family,changed_level)\n'
     t=replace_once(t,first,insert+first)
    else:
     # Boundary signature has an explicit level interval; retain both endpoints.
     pos=t.index('    ',t.index('integer ::')) if False else -1
     lines=t.splitlines(True);i=next(i for i,l in enumerate(lines) if l.lstrip().startswith(('if (','call ')))
     lines.insert(i,'    call op_dependency(3,component,family,first,last)\n');t=''.join(lines)
    s=s[:a]+t+s[b:]
  if s!=original:
   if 'use workflow_profile_195_mod' in s:
    s=re.sub(r'(?m)^  use workflow_profile_195_mod, only:[^\n]*', '  use workflow_profile_195_mod',s)
   else:s=re.sub(r'(?m)^(module \w+[^\n]*\n)',r'\1  use workflow_profile_195_mod\n',s,count=1)
   if 'int64' in s and not re.search(r'use iso_fortran_env[^\n]*int64',s,re.I):
    # Avoid ambiguity if imported transitively by an existing non-private module.
    s=re.sub(r'(?m)^(module \w+[^\n]*\n)',r'\1  use iso_fortran_env, only: int64\n',s,count=1)
   p.write_text(s);changes[p.name]='focused timing/count/dependency hooks; equations unchanged'
 (root/'src/workflow_profile_195.f90').write_bytes(Path(__file__).with_name('workflow_profile_197.f90').read_bytes())
 return changes
