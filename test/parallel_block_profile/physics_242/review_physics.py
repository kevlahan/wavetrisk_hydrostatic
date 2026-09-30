"""Fail closed unless owner-execution qualification is bound to the current source."""
import json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from analyze import analyze
from profile_protocol import checkpoint_name
ROOT=BASE/'.stage242-validation'
tree=ROOT/'candidate'
build=json.loads((tree/'build-242.json').read_text())
assert set(build['binaries'])=={'oracle','production','native','control'}
for name,sha in build['changed_sources'].items():assert digest(BASE/name)==digest(tree/name)==sha,name
for name,sha in build['package_sources'].items():assert digest(tree/name)==digest(BASE/name)==sha,name
for name,sha in build['resolved_sources'].items():assert digest((tree/name).resolve())==sha,name
for name,sha in build['binaries'].items():assert digest(tree/('bin-'+name)/'climate')==sha,name
audit=json.loads((BASE/'docs/parallel-blocks-stage242-dependencies.json').read_text())
for row in audit['sources']:assert digest(BASE/row['file'])==row['sha256']
cases={}
for name,variant,ranks,full in [('oracle-4','oracle',4,False),('production-4','production',4,True),('production-41','production',41,True),('control-4','control',4,False)]:
    out=ROOT/name;r=json.loads((out/'native-state-results.json').read_text());assert r['status']=='passed'
    assert digest(out/'climate')==r['binary_sha256']==build['binaries'][variant]
    for file,key in [('run.log','log_sha256'),('simple.in','input_sha256')]:assert digest(out/file)==r[key]
    rows=r['owner_rows']
    waves=r['batch_waves']
    if variant=='control':assert not r['native_state']
    else:
        states=r['native_state']
        assert sum(c+a for rank,c,a,k in states)==sum(n for rank,n,k in r['physics_columns'])
        assert sum(k for rank,c,a,k in states)==(sum(c+a for rank,c,a,k in states) if variant=='oracle' else 0)
        assert sum(a for rank,c,a,k in states)>0
        for rank in range(ranks):assert sum(row[0]==rank for row in states)==r['steps']
        assert sum(c for rank,c,a in r['native_inputs'])>0
        assert sum(a for rank,c,a in r['native_inputs'])>0
        for rank in range(ranks):
            assert sum(row[0]==rank for row in r['native_inputs'])==r['steps']
            assert sum(row[0]==rank for row in rows)==r['steps']
            assert sum(row[0]==rank for row in waves)==r['steps']
        assert all(w>1 for rank,w in waves)
        assert sum(remote for rank,n,remote,scaffold,checks in rows)>0
        assert sum(n+scaffold for rank,n,remote,scaffold,checks in rows)==sum(n for rank,n,checks in r['physics_columns'])
        assert sum(row[4] for row in rows)==(sum(n for rank,n,checks in r['physics_columns']) if variant=='oracle' else 0)
        if ranks==41:
            assert all(n>0 and remote==n and scaffold==0 for rank,n,remote,scaffold,checks in rows if rank==40)
            assert any(rank==40 and domains==0 and owned>0 for rank,geo,owned,aliases,domains in r['plans'])
        if variant=='oracle':
            assert len(r['exact_alias_checks'])==16
            assert all(checks==30*n for rank,n,checks in r['physics_columns'])
    reference=BASE/('.stage241-validation/production-'+str(ranks))
    observed=[s['state'] for s in analyze(out/'run.log')['step_records']]
    expected=[s['state'] for s in analyze(reference/'run.log')['step_records']]
    assert observed==(expected if full else expected[:1])
    result=dict(status='passed',ranks=ranks,steps=r['steps'],exact_step_states=True,
                owner_columns=sum(row[1] for row in rows),remote_columns=sum(row[2] for row in rows),
                scaffold_columns=sum(row[3] for row in rows),column_oracles=sum(row[4] for row in rows),
                result_sha256=digest(out/'native-state-results.json'),binary_sha256=r['binary_sha256'],
                log_sha256=r['log_sha256'],input_sha256=r['input_sha256'])
    if full:
        f=json.loads((out/'results.json').read_text())
        assert f['steps']==12 and f['coverage']['remaps']==2 and f['numerical']['exact_observed']
        assert f['step_states_exact'] and f['vtk_exact_members']==62
        result.update(remaps=2,vtk_exact_members=62,full_result_sha256=digest(out/'results.json'))
    result['dependency_waves']=sorted({w for rank,w in waves})
    result['native_canonical_requests']=sum(c for rank,c,a in r['native_inputs'])
    result['native_alias_requests']=sum(a for rank,c,a in r['native_inputs'])
    result['retained_canonical_columns']=sum(c for rank,c,a,k in r['native_state'])
    result['retained_alias_columns']=sum(a for rank,c,a,k in r['native_state'])
    result['retained_state_checks']=sum(k for rank,c,a,k in r['native_state'])
    cases[name]=result
for mode in ('route','dependency','version','output','poison','soil','native-route','native-epoch','native-value','entry-poison','retained-owner','retained-alias'):
    out=ROOT/('soil-2' if mode=='soil' else 'reject-'+mode);r=json.loads((out/'results.json').read_text())
    assert r['status']=='passed' and (r['returncode']==0)==(mode in ('soil','poison','entry-poison'))
    assert r['expected_failure']==(mode not in ('soil','poison','entry-poison'))
    assert r['domain_backup']==(mode=='entry-poison')
    assert r['original_source_sha256']==digest(tree/'src/physics/physics_simple.f90')
    for file,key in [('climate','binary_sha256'),('run.log','log_sha256'),('injected.f90','injected_source_sha256'),('simple.in','input_sha256')]:assert digest(out/file)==r[key]
    assert digest(out/checkpoint_name((out/'simple.in').read_text()))==r['seed_sha256']
    assert r['diagnostic'] in (out/'run.log').read_text()
    if mode=='retained-alias':assert 'Inject retained alias corruption' in (out/'run.log').read_text()
    if mode in ('soil','poison','entry-poison'):
        observed=[s['state'] for s in analyze(out/'run.log')['step_records']]
        reference=BASE/'.stage241-validation/production-4'
        n=2 if mode=='soil' else 1
        expected=[s['state'] for s in analyze(reference/'run.log')['step_records']][:n]
        assert len(observed)==n and observed==expected
    cases['reject-'+mode]=r
# Prior frozen numerical evidence is never refreshed to match this change.
prior=BASE/'.stage241-validation/candidate'
old=json.loads((prior/'build-241.json').read_text())
for name,sha in old['resolved_sources'].items():assert digest((prior/name).resolve())==sha,name
for name,sha in old['binaries'].items():assert digest(prior/('bin-'+name)/'climate')==sha,name
main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE,text=True).strip()
assert main=='6a245db7f1632b0205c2bebe95e5f0cd535c2582'
save(BASE/'docs/parallel-blocks-stage242-local-results.json',dict(stage=242,status='validated',cases=cases,
 build=build,dependency_audit_sha256=digest(BASE/'docs/parallel-blocks-stage242-dependencies.json'),
 helpers={str(p.relative_to(BASE)):digest(p) for p in list(Path(__file__).parent.glob('*.py'))+[Path(__file__).with_name('entry_poison.f90.inc'),BASE/'test/parallel_block_profile/physics_240/reject_invalid.py',BASE/'test/parallel_block_profile/physics_239/reject_invalid.py']},
 checkpoint_commit=audit['checkpoint_commit'],main_unchanged=main,prior_frozen_evidence_unchanged=True,
 local_execution=dict(usable_cores=4,solver_jobs_at_once=1,numerical_library_threads=1,rank41_oversubscribed=True),
 default_enabled=False,complete_physics_owner_migration=False,performance_claim=None))
print('PASS Stage242: retained owner physics state, exact 4/41, replay and storage fault checks')
