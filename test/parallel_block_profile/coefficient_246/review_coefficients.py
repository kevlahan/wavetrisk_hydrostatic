"""Bind direct coefficient transfer to its complete reference trajectory."""
import json,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
root=BASE/'.stage246-validation';tree=root/'candidate';out=root/'production-4'
build=json.loads((tree/'build-246.json').read_text())
for n,sha in build['changed_sources'].items():assert digest(BASE/n)==digest(tree/n)==sha,n
for n,sha in build['resolved_sources'].items():assert digest((tree/n).resolve())==sha,n
r=json.loads((out/'lifecycle-results.json').read_text());f=json.loads((out/'results.json').read_text())
assert r['status']=='passed' and r['steps']==12 and r['ranks']==4
assert r['binary_sha256']==build['binaries']['production']==digest(out/'climate')==digest(tree/'bin-production/climate')
assert r['log_sha256']==digest(out/'run.log')
assert f['step_states_exact'] and f['coverage']['remaps']==2 and f['vtk_exact_members']==62 and f['numerical']['passed']
save(BASE/'docs/parallel-blocks-stage246-results.json',dict(stage=246,status='passed',build=build,
    cases={'production-4':dict(path=str((out/'lifecycle-results.json').relative_to(BASE)),sha256=digest(out/'lifecycle-results.json'))},
    scope='Direct vector coefficient transfer, native compression and alias seed. Final independent oracle and 41-rank coverage belong to integrated storage qualification.'))
print('PASS coefficient source-bound review')
