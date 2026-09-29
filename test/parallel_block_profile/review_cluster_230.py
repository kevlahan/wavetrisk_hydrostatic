#!/usr/bin/env python3
"""Recheck collected Stage230 evidence against the sealed campaign and raw logs.

Raw checkpoint/VTK binaries are retained on bbserv; this reviewer independently
checks their collected comparison reports and payload hashes, not omitted data.
"""
import argparse,json
from pathlib import Path
from capability_223 import digest,save,verify_kit
from collect_cluster_logs_230 import verify_archive
from analyze import analyze
from build_capability_230 import source_key
from cluster_230 import campaign
from protocol_230 import ordinary
from validation_230 import storage_evidence
from placement_230 import plan,verify
from rank_metrics_190 import rank_resources


def review(archive,raw,kit,output):
    inventory=verify_archive(archive);verify_kit(kit)
    assert inventory['stage']==230 and inventory['job']=='55089'
    for name,record in inventory['files'].items():
        assert digest(raw/name)==record['sha256'],name
    root=raw/'results';r=json.loads((root/'results.json').read_text())
    manifest=json.loads((kit/'stage230-sources.json').read_text());build=r['build_notes']
    assert r['stage']==230 and r['status']=='passed' and len(r['runs'])==15
    assert build==json.loads((raw/'builds/wt-stage230/capability-230-builds.json').read_text())
    assert build['status']=='passed' and build['manifest_sha256']==digest(kit/'stage230-sources.json')
    assert len(build['physics_objects'])==23
    assert set(build['builds'])=={'J5-baseline','J5-candidate','J5-checked'}
    for key,b in build['builds'].items():
        kind=key.split('-',1)[1]
        assert b['sources']==manifest[source_key(kind)]
        assert digest(raw/f'builds/wt-stage230/build-230-{key}.log')==b['log_sha256']
        for flag in ('FP_CONTRACT=off','STAGE202=on','STAGE203=on','DEBUG='+('check' if kind=='checked' else 'false')):
            assert flag in b['command']
    top=json.loads((root/'topology.json').read_text());cases=plan(top)
    assert top['host']==build['build_host']['node'] and top['host'].split('.')[0]!='bb04'
    for n,case in cases.items():
        assert case==r['cases'][str(n)]
        verify(json.loads((root/f'placement-{n}-verified.json').read_text()),case)
    identity=r['inputs']['J5'];assert digest(root/'inputs/J5/simple.in')==identity['input_sha256']
    rows={};states=None
    for label in r['schedule']:
        row=r['runs'][label];directory=root/label;n=row['ranks'];kind=row['kind'];oracle=row['oracle']
        log=(directory/'run.log').read_text();a=analyze(directory/'run.log')
        recorded=json.loads((directory/'analysis.json').read_text())
        assert Path(recorded['log']).parent.name==label
        assert {k:v for k,v in a.items() if k!='log'}=={k:v for k,v in recorded.items() if k!='log'}
        assert row['returncode']==0 and row['domains']==160 and row['reloads']==[3,4] and row['remaps']==2
        assert row['binary_sha256']==build['builds']['J5-'+kind]['sha256']
        assert digest(directory/'run.log')==row['log_sha256']
        assert digest(directory/'simple.in')==identity['input_sha256']
        resources=rank_resources(directory);verify(resources,cases[n]);assert resources==row['rank_metrics']
        assert len(a['step_records'])==12
        current=[s['state'] for s in a['step_records']]
        if states is None:states=current
        assert current==states
        assert json.loads(json.dumps(storage_evidence(log,n,kind,oracle)))==row['storage']
        if n>160:assert f'Resident block execution: MPI tasks = {n}, geometric Domains = 160' in log
        if oracle:
            for name in ('dynamics','adaptation','vertical-remap'):
                assert f'Block {name} validation oracle enabled' in log
            assert log.count('remap validation oracle max time')==2
        timing=ordinary(log,a,(directory/'test_log').read_text())
        assert timing==row['ordinary'] and timing['median']==row['ordinary_step_median']
        assert [v['index'] for v in timing['excluded']]==[0,7,8]
        reference=root/row['vtk']['reference']
        members=json.loads((directory/'vtk-members.json').read_text())
        assert len(members)==62 and members==json.loads((reference/'vtk-members.json').read_text())
        assert r['runs'][reference.name]['ranks']==n
        assert row['checkpoint_sha256']==r['runs'][reference.name]['checkpoint_sha256']
        for key in (label,'same-count-'+label):
            gate=r['gates'][key];report=gate['report']
            assert gate['passed'] and not gate['failures']
            if key.startswith('same-count-'):
                assert gate['exact_required'] and gate['exact_observed']
                assert report['header_differences']==report['topology_differences']==report['threshold_max_abs']==report['threshold_nonfinite']==0
                assert all(f['different']==f['nonfinite']==f['max_abs']==0 for f in report['fields'].values())
        rows[label]=dict(ranks=n,kind=kind,oracle=oracle,phase=row['phase'],ordinary_median=timing['median'],
                         storage=row['storage'],vtk_exact_members=62,checkpoint_sha256=row['checkpoint_sha256'])
    schedule=[]
    def execute(level,n,kind,oracle,phase,suffix):
        label=f'{level}-n{n}-{phase}-{suffix}-{kind}';schedule.append(label)
        return rows[label]['ordinary_median']
    decisions=campaign(execute,lambda d:None)
    assert schedule==r['schedule'] and decisions==r['decisions']
    confirmation=decisions['confirmation'];assert confirmation['decision']=='within_regression_guard'
    values=[rows[n]['ordinary_median'] for n in schedule if rows[n]['phase']=='confirm']
    result=dict(stage=230,job='55089',archive_sha256=digest(archive),archive_files_verified=len(inventory['files']),
                status='passed',reference_commit=manifest['source_commit'],source_manifest_sha256=digest(kit/'stage230-sources.json'),
                host=top['host'],runs=rows,decisions=decisions,
                confirmation_pair_gains_seconds=[values[0]-values[1],values[3]-values[2]],
                independent_rechecks=['archive and extracted hashes','sealed source/build identities','raw log timing and exclusions',
                    'collective storage census','all-rank placement and resource files','step states and oracle markers',
                    'same-rank VTK member hashes','exact checkpoint comparison reports and hashes','schedule and performance decisions'],
                interpretation='Qualified Domain trend-storage cleanup, within 3% regression guard. Small average timing gain has opposite-sign paired changes; no demonstrated speedup.',
                limitation='Raw checkpoint/VTK files and executables are not in the small review archive; their collected hashes/reports were checked, not recomputed from omitted payloads.',
                next_step='Commit qualified Stage230 cleanup on investigation branch; then migrate full inverse/shared-primitives/tendencies/RK/wavelet producer-consumer cycle to block-owned state. No further trend-storage timing campaign.')
    save(output,result)
    print(json.dumps({k:v for k,v in result.items() if k!='runs'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ('archive','raw','kit','output'):p.add_argument('--'+n,type=Path,required=True)
    review(**vars(p.parse_args()))
