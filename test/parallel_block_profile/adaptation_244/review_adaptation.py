"""Bind native topology-transition qualification to immutable source and binaries."""
import json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
ROOT=BASE/'.stage244-validation';tree=ROOT/'candidate3'
build=json.loads((tree/'build-244.json').read_text())
assert set(build['binaries'])=={'native','production','control','oracle'}
for key in ['changed_sources','package_sources']:
    for n,sha in build[key].items():assert digest(BASE/n)==digest(tree/n)==sha,n
for n,sha in build['resolved_sources'].items():assert digest((tree/n).resolve())==sha,n
for n,sha in build['binaries'].items():assert digest(tree/('bin-'+n)/'climate')==sha,n
cases={}
for name,variant,ranks,checks in [('oracle3-4','native',4,True),('production3-4','production',4,False),('production3-41','production',41,False),('control3-4','control',4,False),('rk-oracle3-4','oracle',4,True)]:
    out=ROOT/name;p=out/'adaptation-results.json';r=json.loads(p.read_text())
    if name=='rk-oracle3-4':
        assert r.get('oracles')
        assert (out/'run.log').read_text().count('Retained RK exact aliases rank=')==16
    assert r['status']=='passed' and r['ranks']==ranks and r['checks']==checks
    assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
    assert digest(out/'run.log')==r['log_sha256']
    assert digest(out/'soil-results.json')==r['soil_results_sha256']
    assert r['steps']==(12 if name.startswith('production') else 1)
    if name.startswith('production'):
        lifecycle=json.loads((out/'results.json').read_text())
        assert lifecycle['steps']==12 and lifecycle['coverage']['remaps']==2
        assert lifecycle['step_states_exact'] and lifecycle['vtk_exact_members']==62 and lifecycle['numerical']['passed']
    cases[name]=dict(path=str(p.relative_to(BASE)),sha256=digest(p),status='passed',steps=r['steps'])
for mode in ['epoch','generation','poison','owner-readback','owner-readback-41']:
    out=ROOT/('reject-'+mode);p=out/'results.json';r=json.loads(p.read_text())
    assert r['status']=='passed' and r['expected_failure']==(mode in ('epoch','generation'))
    source='src/main.f90' if mode=='poison' else 'src/parallel_block_inverse.f90'
    assert r['original_source_sha256']==digest(tree/source)
    assert r['injected_source_sha256']==digest(out/'injected.f90')
    assert r['binary_sha256']==digest(out/'climate') and r['log_sha256']==digest(out/'run.log')
    if mode.startswith('owner-readback'):
        assert (out/'run.log').read_text().count('Native adaptation owner readback exact rank=')==r['ranks']
    if mode=='poison':
        assert (out/'run.log').read_text().count('Native adaptation poison restored rank=')==48
        assert r['lifecycle']['steps']==12 and r['lifecycle']['vtk_exact_members']==62
        assert r['lifecycle']['step_states_exact'] and r['lifecycle']['numerical']['passed']
    cases['reject-'+mode]=dict(path=str(p.relative_to(BASE)),sha256=digest(p),status='passed')
old=BASE/'.stage243b-validation/candidate';prior=json.loads((old/'build-243.json').read_text())
for n,sha in prior['changed_sources'].items():assert digest(old/n)==sha,n
for n,sha in prior['binaries'].items():assert digest(old/('bin-'+n)/'climate')==sha,n
audit=BASE/'docs/parallel-blocks-stage244-kernels.json';kernels=json.loads(audit.read_text())
assert kernels['status']=='passed' and len(kernels['kernels'])==15
for n,sha in kernels['sources'].items():assert digest(BASE/n)==sha
main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
helpers={str(p.relative_to(BASE)):digest(p) for p in sorted(Path(__file__).parent.iterdir()) if p.is_file()}
save(BASE/'docs/parallel-blocks-stage244-results.json',dict(stage=244,status='passed',main_unchanged=main,build=build,helpers=helpers,cases=cases,kernel_audit_sha256=digest(audit),scope='Native numerical adaptation and topology transition; lifecycle adapters and persistent storage cleanup remain'))
print('PASS Stage244 source-bound reviewer: ten cases',flush=True)
