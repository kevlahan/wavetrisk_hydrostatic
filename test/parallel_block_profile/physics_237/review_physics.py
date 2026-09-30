"""Fail closed unless the bounded physics migration and source identities pass."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze

def review(root):
    root=root.resolve();tree=root/'candidate';build=json.loads((tree/'build-237.json').read_text())
    assert set(build['binaries'])=={'production','native','oracle','control'}
    for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
    for name,sha in build['sources'].items():assert digest(tree/name)==sha,name
    for name,sha in build['resolved_sources'].items():
        assert digest((tree/name).resolve())==sha,name
        if name in build['changed_sources']:assert digest((BASE/name).resolve())==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
    audit=json.loads((BASE/'docs/parallel-blocks-stage237-dependencies.json').read_text())
    for row in audit['dependencies']:assert digest(BASE/row['file'])==row['sha256']
    cases={};reference=BASE/'.stage236-validation/qualified-production-4'
    first=analyze(reference/'run.log')['step_records'][0]['state']
    for name,variant,full in [('production-4','production',True),('production-41','production',True),('oracle-4','oracle',False),('control-4','control',False)]:
        out=root/name;r=json.loads((out/'physics-results.json').read_text());assert r['status']=='passed'
        assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
        for file,key in [('run.log','log_sha256'),('simple.in','input_sha256')]:assert digest(out/file)==r[key]
        if full:
            f=json.loads((out/'results.json').read_text())
            assert f['steps']==12 and f['coverage']['remaps']==2 and f['numerical']['exact_observed']
            assert f['step_states_exact'] and f['vtk_exact_members']==62;r['full_validation']=f
        else:assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[first]
        if variant=='oracle':
            assert len(r['exact_alias_checks'])==16
            assert sum(n for rank,n,checks in r['physics_columns'])>0
            assert all(checks==30*n for rank,n,checks in r['physics_columns'])
        if r['ranks']>40:
            assert any(rank==40 and domains==0 and owned>0 for rank,geo,owned,aliases,domains in r['plans'])
            assert sum(rank==40 and n==0 for rank,n,checks in r['physics_columns'])==12
        for key in ('counts','primitive_counts','pentagon_counts','gradient_counts','rk_closure','exact_alias_checks','physics_columns'):
            r[key+'_records']=len(r[key]);r[key]=list({tuple(row[:2]) if key in ('rk_closure','exact_alias_checks') else row[0]:row for row in r[key]}.values())
        cases[name]=r
    for mode in ('poison','epoch','duplicate','source'):
        out=root/('poison-4' if mode=='poison' else 'reject-'+mode);r=json.loads((out/'results.json').read_text())
        assert r['status']=='passed' and (r['returncode']==0)==(mode=='poison')
        assert r['original_source_sha256']==digest(tree/'src/physics/physics_simple.f90')
        for file,key in [('climate','binary_sha256'),('run.log','log_sha256'),('injected.f90','injected_source_sha256')]:assert digest(out/file)==r[key]
        assert r['diagnostic'] in (out/'run.log').read_text()
        if mode=='poison':
            assert not r['domain_backup']
            assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[first]
        cases[mode]=r
    main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip();assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    save(BASE/'docs/parallel-blocks-stage237-local-results.json',dict(stage=237,status='physics-column-scratch-migration-validated',
         complete_physics_owner_migration=False,complete_domain_storage_removal=False,default_enabled=False,
         checkpoint_commit=audit['checkpoint_commit'],build=build,cases=cases,raw_root=str(root),performance_claim=None,
         local_execution=dict(usable_cores=4,numerical_library_threads=1,rank41_oversubscribed=True,
                              scheduling='Early cases overlapped; remaining cases serialized after user clarification. Correctness evidence only.'),
         dependency_audit_sha256=digest(BASE/'docs/parallel-blocks-stage237-dependencies.json'),
         helpers={str(p.relative_to(BASE)):digest(p) for p in Path(__file__).parent.glob('*.py')},main_unchanged=main))
    print('PASS Stage237 physics scratch: exact 4/41 trajectories, layer oracles, poison and faults')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=BASE/'.stage237-validation');review(p.parse_args().root)
