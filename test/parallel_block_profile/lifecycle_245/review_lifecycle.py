"""Bind between-step lifecycle checks to the frozen candidate."""
import json,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
root=BASE/'.stage245-validation';tree=root/'candidate2'
build=json.loads((tree/'build-245.json').read_text())
for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
for name,sha in build['resolved_sources'].items():assert digest((tree/name).resolve())==sha,name
cases={}
for name,variant,steps in [('production2-4','production',12),('check2-4','native',1)]:
    out=root/name;p=out/'lifecycle-results.json';r=json.loads(p.read_text())
    assert r['status']=='passed' and r['steps']==steps and r['ranks']==4
    assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]==digest(tree/('bin-'+variant)/'climate')
    assert digest(out/'run.log')==r['log_sha256']
    assert digest(out/'adaptation-results.json')==r['adaptation_results_sha256']
    if steps==12:
        full=json.loads((out/'results.json').read_text())
        assert full['step_states_exact'] and full['coverage']['remaps']==2 and full['vtk_exact_members']==62
        assert full['numerical']['passed']
    cases[name]=dict(path=str(p.relative_to(BASE)),sha256=digest(p))
save(BASE/'docs/parallel-blocks-stage245-results.json',dict(stage=245,status='passed',build=build,cases=cases,
    scope='Native between-step lifecycle; legacy coefficient transport and persistent views still require cleanup. Final 41-rank qualification deferred to integrated storage candidate.'))
print('PASS lifecycle source-bound review')
