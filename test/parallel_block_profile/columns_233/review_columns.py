"""Verify the frozen build and local evidence for the partial owner migration."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze

def review(root):
    root=root.resolve();tree=root/'final-candidate'
    build=json.loads((tree/'build-233.json').read_text())
    for name,sha in build['changed_sources'].items():
        assert digest(BASE/name)==sha and digest(tree/name)==sha,name
    for name,sha in build['sources'].items():assert digest(tree/name)==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/f'bin-233-{name}/climate')==sha,name
    cases={}
    for name,variant in [('qualified-oracle-4','oracle'),('final-production-41','production'),('control-4','control')]:
        out=root/name;result=json.loads((out/'columns-results.json').read_text())
        assert result['status']=='passed'
        assert digest(out/'climate')==result['binary_sha256']==build['binaries'][variant],name
        assert digest(out/'run.log')==result['log_sha256']
        assert digest(out/'simple.in')==result['input_sha256']
        if variant=='control':
            assert result['control'] and result['steps']==1 and not result['counts']
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
                assert all(any(r==rank and executed>0 for r,executed,checked in result['counts']) for rank in zero_owners)
        cases[name]=result
    for mode in ('missing','corrupt','epoch','fallback'):
        out=root/('reject-'+mode);result=json.loads((out/'results.json').read_text())
        assert result['status']=='passed' and result['expected_failure'] and result['returncode']!=0
        assert digest(out/'climate')==result['binary_sha256']
        assert digest(out/'run.log')==result['log_sha256']
        assert digest(out/'injected.f90')==result['injected_source_sha256']
        assert result['diagnostic'] in (out/'run.log').read_text()
        cases['reject-'+mode]=result
    assert subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    save(BASE/'docs/parallel-blocks-stage233-local-results.json',dict(stage=233,status='partial-owner-integration-validated',
        complete_tendency_migration=False,default_enabled=False,build=build,cases=cases,
        scope='RK interstage hydrostatic columns execute on final block owners; reference aliases and compatibility result consumers remain.',
        helpers={str(p.relative_to(BASE)):digest(p) for p in Path(__file__).parent.glob('*.py')},
        raw_root=str(root),performance_claim=None,
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=BASE,text=True).strip(),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
        main_unchanged='6a245db7f1632b0205c2bebe95e5f0cd535c2582'))
    print('PASS: two full trajectories, feature-off control and four negative guards; full tendency migration remains incomplete')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=BASE/'.stage233-validation')
    review(p.parse_args().root)
