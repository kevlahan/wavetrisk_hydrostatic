"""Require source-bound RK trajectory, storage-absence, oracle and fault evidence."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze

def review(root):
    root=root.resolve();tree=root/'qualified-candidate'
    build=json.loads((tree/'build-236.json').read_text())
    assert set(build['binaries'])=={'production','native','oracle','control'}
    for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
    for name,sha in build['sources'].items():assert digest(tree/name)==sha,name
    for name,sha in build['resolved_sources'].items():
        assert digest((tree/name).resolve())==sha,name
        if name in build['changed_sources']:assert digest((BASE/name).resolve())==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
    cases={}
    for name,variant,full in [('qualified-production-4','production',True),('qualified-production-41','production',True),
                             ('qualified-all-oracles-4','oracle',False),('qualified-native-4','native',False),
                             ('control-4','control',False)]:
        out=root/name;r=json.loads((out/'rk-results.json').read_text())
        assert r['status']=='passed'
        assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
        for file,key in [('run.log','log_sha256'),('simple.in','input_sha256')]:assert digest(out/file)==r[key]
        if full:
            f=json.loads((out/'results.json').read_text())
            assert f['steps']==12 and f['coverage']['remaps']==2 and f['numerical']['exact_observed']
            assert f['step_states_exact'] and f['vtk_exact_members']==62
            assert f['binary_sha256']==r['binary_sha256'];r['full_validation']=f
        else:
            assert r['steps']==1
            states=[s['state'] for s in analyze(out/'run.log')['step_records']]
            expected=[s['state'] for s in analyze(root/'qualified-production-4/run.log')['step_records']]
            assert states==expected[:1]
        log=(out/'run.log').read_text()
        if variant=='oracle':
            assert len(r['exact_alias_checks'])==16
            for group in ('dynamics','adaptation','vertical-remap'):assert f'Block {group} validation oracle enabled' in log
        elif variant!='control':
            assert len(r['q1_absence'])==r['ranks']*r['steps']
            assert len(r['return_storage_absence'])==4*r['ranks']*r['steps']
        if r['ranks']>40:
            zero={rank for rank,geo,owned,aliases,domains in r['plans'] if domains==0 and owned>0}
            assert 40 in zero
            assert all(any(row[0]==rank and row[3]>0 for row in r['rk_closure']) for rank in zero)
        for key in ('counts','primitive_counts','pentagon_counts','gradient_counts','rk_closure','exact_alias_checks'):
            r[key+'_records']=len(r[key]);r[key]=list({tuple(row[:2]) if key in ('rk_closure','exact_alias_checks') else row[0]:row for row in r[key]}.values())
        cases[name]=r
    for mode in ('missing','epoch','source','base','fallback','poison'):
        out=root/('poison-4' if mode=='poison' else 'reject-'+mode)
        r=json.loads((out/'results.json').read_text())
        original={'base':'parallel_block.f90','fallback':'time_integr.f90'}.get(mode,'parallel_block_inverse.f90')
        assert r['original_source_sha256']==digest(tree/'src'/original)
        assert r['status']=='passed' and r['expected_failure']==(mode!='poison')
        assert (r['returncode']==0)==(mode=='poison')
        for file,key in [('climate','binary_sha256'),('run.log','log_sha256'),('injected.f90','injected_source_sha256')]:
            assert digest(out/file)==r[key]
        assert r['diagnostic'] in (out/'run.log').read_text()
        if mode=='poison':
            assert r['primary_fields_poisoned'] and not r['domain_backup']
            assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[analyze(root/'qualified-production-4/run.log')['step_records'][0]['state']]
        cases[mode]=r
    main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
    assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    save(BASE/'docs/parallel-blocks-stage236-local-results.json',dict(stage=236,status='rk-migration-validated',
        complete_rk_migration=True,complete_domain_storage_removal=False,default_enabled=False,
        configuration='Restarted climate Simple physics RK4',build=build,cases=cases,raw_root=str(root),performance_claim=None,
        helpers={str(p.relative_to(BASE)):digest(p) for p in Path(__file__).parent.glob('*.py')},main_unchanged=main))
    print('PASS Stage236: retained RK, full 4/41 trajectories, independent oracles, poison and producer faults')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=BASE/'.stage236-validation')
    review(p.parse_args().root)
