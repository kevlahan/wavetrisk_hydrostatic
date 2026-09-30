"""Bind accepted soil runs, adversarial checks, helpers and sources to one frozen candidate."""
import json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
ROOT=BASE/'.stage243b-validation';tree=ROOT/'candidate'
build=json.loads((tree/'build-243.json').read_text())
assert set(build['binaries'])=={'native','production','control'}
for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
for name,sha in build['package_sources'].items():assert digest(tree/name)==digest(BASE/name)==sha,name
for name,sha in build['resolved_sources'].items():assert digest((tree/name).resolve())==sha,name
for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
cases={}
for name,variant,ranks,check in [('oracle-4','native',4,True),('production-4','production',4,False),('production-41','production',41,False),('control-4','control',4,False)]:
    out=ROOT/name;r=json.loads((out/'soil-results.json').read_text())
    assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
    assert digest(out/'run.log')==r['log_sha256']
    assert r['soil_check']==check
    assert r['steps']==(12 if name.startswith('production') else 1)
    if name.startswith('production'):
        lifecycle=json.loads((out/'results.json').read_text())
        assert lifecycle['steps']==12 and lifecycle['coverage']['remaps']==2
        assert lifecycle['step_states_exact'] and lifecycle['vtk_exact_members']==62
    cases[name]=r
for mode in ['epoch','range','alias','poison','owner-readback','owner-readback-41']:
    out=ROOT/('reject-poison-exact-capacity' if mode=='poison' else 'reject-'+mode)
    r=json.loads((out/'results.json').read_text())
    assert r['status']=='passed' and r['expected_failure']==(mode not in ('poison','owner-readback','owner-readback-41'))
    assert r['original_source_sha256']==digest(tree/'src/parallel_block_inverse.f90')
    assert r['injected_source_sha256']==digest(out/'injected.f90')
    assert r['binary_sha256']==digest(out/'climate') and r['log_sha256']==digest(out/'run.log')
    if mode.startswith('owner-readback'):
        assert (out/'run.log').read_text().count('Native soil owner readback exact rank=')==r['ranks']
    if mode=='poison':
        assert (out/'run.log').read_text().count('Native soil Domain poison restored after transform rank=')==48
        assert r['lifecycle']['steps']==12 and r['lifecycle']['vtk_exact_members']==62
        assert r['lifecycle']['step_states_exact'] and r['lifecycle']['numerical']['passed']
    cases['reject-'+mode]=r
old=BASE/'.stage242-validation/candidate';prior=json.loads((old/'build-242.json').read_text())
for n,sha in prior['changed_sources'].items():assert digest(old/n)==sha,n
for n,sha in prior['binaries'].items():assert digest(old/('bin-'+n)/'climate')==sha,n
audit=BASE/'docs/parallel-blocks-stage243-kernels.json'
kernels=json.loads(audit.read_text())
assert kernels['status']=='passed' and len(kernels['kernels'])==15
for n,sha in kernels['sources'].items():assert digest(BASE/n)==sha
main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
helpers={str(p.relative_to(BASE)):digest(p) for p in sorted(Path(__file__).parent.iterdir()) if p.is_file()}
save(BASE/'docs/parallel-blocks-stage243-results.json',dict(stage=243,status='passed',main_unchanged=main,build=build,helpers=helpers,cases=cases,kernel_audit_sha256=digest(audit),
    scope='Native soil transform/inverse with exact alias traversal and selected-range owner storage; adaptation and bounded lifecycle adapters remain'))
print('PASS Stage243 source-bound reviewer: ten cases',flush=True)
