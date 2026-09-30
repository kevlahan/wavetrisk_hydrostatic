#!/usr/bin/env python3
"""Source-anchored migration ledger; classifications are reviewed, not inferred liveness."""
import hashlib,json,re,subprocess
from pathlib import Path
from capability_223 import save
BASE=Path(__file__).resolve().parents[2]
ENTRIES=[
 ('surface/column primitives','src/diagnostics.f90',['cal_surf_press','integrate_pressure_up'],'Domain owner','q/sol_mean -> surf_press, press, geopot, exner','Move ordered full-column producer and its actual primitive readers together; existing cached block Exner is not the restricted primitive Exner.'),
 ('shared primitive and diffusion inputs','src/multi_level.f90',['block_tendency_compatibility_ml','basic_operators','cal_Laplacian_scalars'],'Domain owner','q, sol_mean, Domain geometry and diffusion/primitive scratch','Compile block-local operand/alias topology; preserve descending horizontal levels and ascending physical layers.'),
 ('mass/velocity producer','src/multi_level.f90',['compute_native_velocity_source','compute_native_velocity_gradient'],'Domain owner','native_mass(d), native_velocity(d), Domain PV/Bernoulli/Exner and physics callbacks','Migrate stencil, restriction and dependency exchange with producer output; no more geometry-indexed replacement arrays.'),
 ('retained native velocity delivery','src/parallel_block_mpi.f90',['prepare_block_native_velocity_remainder','refresh_candidate_block_velocity_remainder'],'Domain owner -> block owner','native velocity output; production Domain solution argument now absent','Delete source-owner return after its producer executes on final block owners. Current transport is still real.'),
 ('retained tendency assembly','src/parallel_block_mpi.f90',['capture_block_domain_multistage_candidate_tendency','evaluate_candidate_block_velocity_recomposition','compute_block_exner_difference_kernel'],'block owner','retained scalar and velocity tendencies; compact scalar columns; production Domain solution argument now absent','Keep true inputs and stage-base checkpoint. Preserve add/subtract rounding contract until explicitly validated change.'),
 ('RK and wavelet stage','src/time_integr.f90',['RK4','RK_sub_step_compatibility','native_provisional_WT'],'mixed Domain/block owners','Domain q1/base q plus native geometry-indexed RHS and retained block RK state','Replace complete stage input/output path, including untouched scaffold, before deleting q1 or sol.'),
 ('inverse aliases','src/parallel_block_inverse.f90',['publish_native_inverse'],'geometry route owner -> Domain consumers','canonical workspace values and phase-specific aliases','Preserve storage/version identity for scaffold/poles/aliases, not just coincident geometric coordinates.'),
 ('external physics and output','src/main.f90',['apply_domain_physics_consumers','update_domain_soil_wavelet_consumer','prepare_end_step_domain_consumers'],'Domain owner','Simple physics, soil, adaptation boundary and output/checkpoint image','Migrate supported physics/soil state and retain only bounded explicit external adapters. Geometry/connectivity are not prognostic storage.'),
 ('unused derived hydrostatic cache','src/parallel_block.f90',['compute_local_block_hydrostatic_patch','apply_local_block_hydrostatic_consumer','local_block_hydrostatic_statistics'],'block owner, validation only','surface_pressure, dynamic_exner, air_temperature','Removed from production readiness/allocation/refresh; explicit validation retains it. Not the live primitive pressure calculation.'),
]

def audit():
 rows=[]
 for family,file,names,owner,state,obligation in ENTRIES:
  text=(BASE/file).read_text();anchors=[]
  for name in names:
   m=re.search(r'^\s*(?:recursive\s+)?subroutine\s+'+name+r'\b',text,re.M|re.I)
   if m is None:raise ValueError((file,name))
   end=re.search(r'^[ \t]*end subroutine(?:[ \t]+'+name+r')?[ \t]*$',text[m.end():],re.M|re.I)
   if end is None:raise ValueError(('end',file,name))
   body=text[m.start():m.end()+end.end()]
   anchors.append(dict(routine=name,line=text[:m.start()].count('\n')+1,sha256=hashlib.sha256(body.encode()).hexdigest()))
  rows.append(dict(family=family,file=file,execution_owner=owner,state=state,required_migration=obligation,routines=anchors))
 result=dict(stage=231,base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
             classification='Manual producer/consumer review with executable source-location/hash verification; not dynamic timing attribution.',
             dependencies=rows,
             scaling={'source_fact':'Loops over size(grid) execute no Domain primitive/RK work on zero-Domain ranks. Native mass/velocity workspaces retain this ownership. Block count is determined by available four-child subtrees, not by MPI task count.',
                      'inference':'Remaining owner-limited work and Domain/block transfers can limit scaling past the Domain count; the relative contribution of computation, waiting, balance and bandwidth is not measured by this source audit.',
                      'follow_up':'After migration, compare fixed workload/block topology and placement at task counts below/above Domain count. Use bounded per-rank phase compute/wait evidence and separate-run confirmation, not timestep pseudoreplication.'})
 save(BASE/'docs/parallel-blocks-stage231-dependencies.json',result)
 print('PASS:',len(rows),'reviewed dependency families anchored to',sum(len(r['routines']) for r in rows),'source routines')
if __name__=='__main__':audit()
