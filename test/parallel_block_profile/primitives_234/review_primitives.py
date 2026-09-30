"""Verify the frozen build and local evidence for the partial owner migration."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze

def review(root):
    root=root.resolve();tree=root/'qualified-candidate'
    build=json.loads((tree/'build-234.json').read_text())
    for name,sha in build['changed_sources'].items():
        assert digest(BASE/name)==sha and digest(tree/name)==sha,name
    for name,sha in build['sources'].items():assert digest(tree/name)==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/f'bin-234-{name}/climate')==sha,name
    cases={}
    for name,variant in [('qualified-oracle-4','oracle'),('qualified-production-41','production'),('control-4','control')]:
        out=root/name;result=json.loads((out/'primitives-results.json').read_text())
        assert result['status']=='passed'
        assert digest(out/'climate')==result['binary_sha256']==build['binaries'][variant],name
        assert digest(out/'run.log')==result['log_sha256']
        assert digest(out/'simple.in')==result['input_sha256']
        if variant=='control':
            assert result['control'] and result['steps']==1 and not result['counts'] and not result['retention']
            states=[s['state'] for s in analyze(out/'run.log')['step_records']]
            expected=[s['state'] for s in analyze(root/'qualified-oracle-4/run.log')['step_records']]
            assert states==expected[:1]
        else:
            full=json.loads((out/'results.json').read_text())
            assert full['steps']==12 and full['coverage']['remaps']==2
            assert full['numerical']['exact_observed'] and full['step_states_exact'] and full['vtk_exact_members']==62
            assert full['binary_sha256']==result['binary_sha256']
            result['full_validation']=full
            if variant=='oracle':assert all(checked>0 for rank,executed,checked in result['counts'])
            if variant=='production':
                assert all(checked==0 for rank,executed,checked in result['counts'])
                zero_owners={rank for rank,geo,owned,aliases,domains in result['plans'] if owned>0 and domains==0}
                assert zero_owners
                assert all(any(r==rank and executed>0 for r,executed,checked in result['gradient_counts']) for rank in zero_owners)
                assert all(any(r==rank and executed>0 for r,executed,checked in result['counts']) for rank in zero_owners)
        cases[name]=result
    out=root/'qualified-all-oracles-4'
    result=json.loads((out/'primitives-results.json').read_text())
    assert result['status']=='passed' and result['steps']==1
    assert digest(out/'climate')==result['binary_sha256']==build['binaries']['oracle']
    assert digest(out/'run.log')==result['log_sha256']
    log=(out/'run.log').read_text()
    for group in ('dynamics','adaptation','vertical-remap'):
        assert f'Block {group} validation oracle enabled' in log
    assert 'Domain tendency storage: legacy/oracle' in log
    expected=analyze(root/'qualified-oracle-4/run.log')['step_records'][0]['state']
    assert analyze(out/'run.log')['step_records'][0]['state']==expected
    cases[out.name]=result
    for name,result in cases.items():
        if result['control']:continue
        for family,width in [('primitive_counts',2),('pentagon_counts',1),('gradient_counts',1)]:
            last={r[0]:r[1:] for r in result[family]}
            total=[sum(row[j] for row in last.values()) for j in range(2*width)]
            if 'production' in name:
                assert total[width:]==[0]*width,(name,family,total)
            else:
                assert total[:width]==total[width:] and all(x>0 for x in total[:width]),(name,family,total)
        assert {r for r,n in result['retention']}==set(range(result['ranks']))
    for mode in ('primitive','pentagon','gradient','missing','epoch','capture'):
        out=root/('reject-'+mode);result=json.loads((out/'results.json').read_text())
        assert result['status']=='passed' and result['expected_failure'] and result['returncode']!=0
        assert digest(out/'climate')==result['binary_sha256']
        assert digest(out/'run.log')==result['log_sha256']
        assert digest(out/'injected.f90')==result['injected_source_sha256']
        assert result['diagnostic'] in (out/'run.log').read_text()
        cases['reject-'+mode]=result
    assert subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    # Keep the durable report compact; the raw log hash binds every comparison.
    for result in cases.values():
        if result.get('expected_failure'):continue
        for key in ('counts','primitive_counts','pentagon_counts','gradient_counts'):
            rows=result[key]
            result[key+'_records']=len(rows)
            result[key]=list({row[0]:row for row in rows}.values())
        rows=result['retention']
        result['retention_records']=len(rows)
        result['retention']=list({row[0]:row for row in rows}.values())
    save(BASE/'docs/parallel-blocks-stage234-local-results.json',dict(stage=234,status='partial-owner-integration-validated',
        complete_tendency_migration=False,default_enabled=False,build=build,cases=cases,
        scope='RK interstage regular/SW/pentagon primitives and final gradients execute on block owners, with direct retained velocity installation; full tendency migration remains incomplete.',
        helpers={str(p.relative_to(BASE)):digest(p) for p in Path(__file__).parent.glob('*.py')},
        raw_root=str(root),performance_claim=None,
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=BASE,text=True).strip(),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
        main_unchanged='6a245db7f1632b0205c2bebe95e5f0cd535c2582'))
    print('PASS: two full trajectories, complete one-step oracles, feature-off control and six injected faults; full tendency migration remains incomplete')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=BASE/'.stage234-validation')
    review(p.parse_args().root)
