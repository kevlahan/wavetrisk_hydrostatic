"""Source-anchored read/write contract for the Simple-physics velocity boundary."""
import re,sys,subprocess
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
SOURCES={
 'src/diagnostics.f90':['kinetic_energy'],
 'src/domain.f90':['id_edge'],
 'src/domain_ops.f90':['apply_no_bdry2'],
 'src/physics/physics_simple.f90':['physics_simple_step','plan_velocity_column','physics_call'],
 'src/parallel_block_primitive_kernels.f90':['physics_kinetic_energy'],
 'src/main.f90':['apply_domain_physics_consumers','update_domain_soil_wavelet_consumer','time_step']}
rows=[]
for name,routines in SOURCES.items():
    text=(BASE/name).read_text();anchors=[]
    for routine in routines:
        match=re.search(r'^.*(?:subroutine|function)\s+'+routine+r'\b',text,re.M|re.I)
        assert match,(name,routine)
        anchors.append(dict(routine=routine,line=text[:match.start()].count('\n')+1))
    rows.append(dict(file=name,sha256=digest(BASE/name),anchors=anchors))
save(BASE/'docs/parallel-blocks-stage238-dependencies.json',dict(stage=238,
 checkpoint_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),sources=rows,
 operands=[dict(slot=n+1,offset=offset,edge=edge) for n,(offset,edge) in enumerate([('current','UP'),('current','DG'),('current','RT'),('south','UP'),('southwest','DG'),('west','RT')])],
 read_version='Zero: local alias at epoch entry. Positive: ordinal of latest preceding non-pole writer in the existing walk.',
 writes='Non-pole columns write all three id_edge velocities at every atmospheric layer; poles do not write velocity. Temperature and soil updates remain in place.',
 ordering='apply_no_bdry2: ascending levels, ordinary columns then poles at each level. Exact full local order retained, including read-before-write hazards.',
 aliases='Storage identities are Domain-local. An alias with no local writer retains version zero even if another copy has a writer. Do not replace it with a newer canonical-owner value.',
 epochs='Prepare after final RK Domain publication; destroy before soil wavelets/adaptation. Pending atmospheric wavelet transaction and soil lifetime are unchanged.',
 next_boundary='Map these local storage-version pairs to canonical owner operands without refreshing entry aliases; assess package mutable state before dispatch.',
 complete_physics_owner_migration=False,performance_claim=None))
print('PASS velocity dependency audit')
