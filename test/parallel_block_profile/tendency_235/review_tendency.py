"""Bind Stage235 source, binaries and correctness evidence before reporting completion."""
import argparse,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze

def review(root):
    root=root.resolve();tree=root/'qualified-candidate'
    build=json.loads((tree/'build-235.json').read_text())
    for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
    for name,sha in build['sources'].items():assert digest(tree/name)==sha,name
    for name,sha in build['source_contents'].items():
        assert digest((tree/name).resolve())==sha,name
        if name in build['changed_sources']:assert digest((BASE/name).resolve())==sha,name
    for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
    cases={}
    for name,variant,full in [('qualified-production-4','production',True),('qualified-production-41','production',True),
                             ('qualified-all-oracles-4','oracle',False),
                             ('control-4','control',False)]:
        out=root/name;r=json.loads((out/'tendency-results.json').read_text())
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
        if variant=='oracle':
            assert r['source_checks'] and sum(n for rank,n in r['source_checks'])>0
            assert r['mass_checks'] and sum(n for rank,n in r['mass_checks'])>0
        if r['ranks']>40:
            zero={rank for rank,geo,owned,aliases,domains in r['plans'] if domains==0 and owned>0}
            assert zero
            assert all(any(row[0]==rank and row[2]>0 for row in r['owner_programs']) for rank in zero)
        if name=='qualified-all-oracles-4':
            log=(out/'run.log').read_text()
            for group in ('dynamics','adaptation','vertical-remap'):assert f'Block {group} validation oracle enabled' in log
        # Preserve full logs by hash; keep only final per-rank counters in report.
        for key in ('counts','primitive_counts','pentagon_counts','gradient_counts','mass_checks','source_checks'):
            r[key+'_records']=len(r[key]);r[key]=list({row[0]:row for row in r[key]}.values())
        cases[name]=r
    for mode in ('missing','epoch','source','scratch'):
        out=root/('scratch-4' if mode=='scratch' else 'reject-'+mode)
        r=json.loads((out/'results.json').read_text())
        original={'missing':'parallel_block_diffusion.f90','epoch':'parallel_block_velocity_owner.f90',
                  'source':'parallel_block_program.f90','scratch':'multi_level.f90'}[mode]
        assert r['original_source_sha256']==digest(tree/'src'/original)
        assert r['status']=='passed' and r['expected_failure']==(mode!='scratch')
        assert (r['returncode']==0)==(mode=='scratch')
        for file,key in [('climate','binary_sha256'),('run.log','log_sha256'),('injected.f90','injected_source_sha256')]:
            assert digest(out/file)==r[key]
        assert r['diagnostic'] in (out/'run.log').read_text()
        if mode=='scratch':
            assert r['primary_fields_poisoned'] and r['removed_scratch_poisoned']
            assert [s['state'] for s in analyze(out/'run.log')['step_records']]==[analyze(root/'qualified-production-4/run.log')['step_records'][0]['state']]
        cases[mode]=r
    main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
    assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
    save(BASE/'docs/parallel-blocks-stage235-local-results.json',dict(stage=235,status='all-stage-tendency-kernels-validated',
        complete_tendency_kernel_migration=True,complete_domain_storage_removal=False,default_enabled=False,
        build=build,cases=cases,raw_root=str(root),performance_claim=None,
        helpers={str(p.relative_to(BASE)):digest(p) for p in Path(__file__).parent.glob('*.py')},main_unchanged=main))
    print('PASS Stage235: all-stage kernels, full 4/41 trajectories, independent oracles, scratch poison and faults')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=BASE/'.stage235-validation')
    review(p.parse_args().root)
