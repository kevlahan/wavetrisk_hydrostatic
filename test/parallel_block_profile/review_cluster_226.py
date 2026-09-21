import sys,json,re,hashlib,subprocess
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'test/parallel_block_profile'))
from capability_226 import events,schedule,digest
from analyze import analyze
from checkpoint_gate import accept
from placement_226 import plan,verify
from rank_metrics_190 import rank_resources
from collect_cluster_logs_226 import verify_archive
base=Path('/private/tmp/stage226-review-54731');root=base/'results'
archive=Path('/Users/kevlahan/Desktop/stage226-review-54731.tar.gz')
inventory=verify_archive(archive)
r=json.loads((root/'results.json').read_text())
kit=Path('/private/tmp/wavetrisk-stage226-cluster-kit')
m=json.loads((kit/'stage226-sources.json').read_text())
for name in ('src/comm_mpi.f90','test/climate/climate.f90','test/climate/test_case_module.f90'):
    assert digest(Path(name))==m['candidate_sources'][name],name
assert r['stage']==226 and r['status']=='passed' and inventory['job']=='54731'
expected=[f'{l}-n{n}-{k}-'+('oracles' if o else 'production')+'-'+lay for l,n,k,o,lay in schedule()]
assert r['schedule']==expected and set(r['runs'])==set(expected)
b=r['build_notes'];assert b['status']=='passed' and b['manifest_sha256']==digest(kit/'stage226-sources.json')
assert b['fp_contract']=='off' and len(b['builds'])==6 and len(b['physics_objects'])==23
for key,v in b['builds'].items():
    assert v['sources']==m['reference_sources' if key.endswith('baseline') else key.split('-')[1]+'_sources']
    assert v['log_sha256']==digest(base/'builds/wt-stage226'/f'build-226-{key}.log')
    assert 'DEBUG='+('check' if key.endswith('checked') else 'false') in v['command']
    assert all(flag in v['command'] for flag in ['FP_CONTRACT=off','STAGE202=on','STAGE203=off'])
top=json.loads((root/'topology.json').read_text());plans=plan(top)
assert b['build_host']['node']==top['host']
for n,p in plans.items():
    assert p==r['cases'][str(n)]
    verify(json.loads((root/f'placement-{n}-verified.json').read_text()),p)
rows={};trajectories={};checksums={}
for label,row in r['runs'].items():
    d=root/label;log=(d/'run.log').read_text();a=analyze(d/'run.log')
    assert row['returncode']==0 and a['completed'] and len(a['step_records'])==row['steps']==12
    assert digest(d/'run.log')==row['log_sha256']
    key=row['level']+'-'+row['kind'];assert row['binary_sha256']==b['builds'][key]['sha256']
    verify(rank_resources(d),plans[row['ranks']])
    reloads=list(map(int,re.findall(r'Restarting from checkpoint\s+(\d+)',log)))
    assert log.count('Saving field')==2 and log.count('Saving checkpoint')==1
    assert reloads==row['reloads']==([2,3] if row['level']=='J4' else [3,4])
    assert log.count('Remapping vertical coordinates')==row['remaps']==2
    trajectory=[s['state'] for s in a['step_records']]
    level=row['level'];trajectories.setdefault(level,trajectory);assert trajectory==trajectories[level]
    checksums.setdefault(level,row['checkpoint_sha256']);assert row['checkpoint_sha256']==checksums[level]
    summary={k:row[k] for k in ['ranks','level','kind','oracle','layout','steps','reloads','remaps','checkpoint_sha256']}
    if row['kind']=='checked':
        ev=events(d,row['ranks'],row['domains'],12,row['layout']);assert ev==row['capability']
        zero=sorted({tuple(z) for z in ev['zero_domain_rank_sets']})
        summary.update(zero_domain_rank_sets=zero,restart_source_ranks=ev['restart_source_ranks'])
        emptycounts=[]
        for rank in set(sum((list(z) for z in zero),[])):
            for line in (d/f'capability223-rank-{rank}.txt').read_text().splitlines():
                f=line.split()
                if f[0]!='restart' and int(f[3])==0:emptycounts.append(int(f[4]))
        if emptycounts:summary['zero_domain_owned_blocks_range']=[min(emptycounts),max(emptycounts)]
    if row['kind']=='candidate':
        assert not list(d.glob('capability223-rank-*.txt'))
        assert not any(k.startswith('WAVETRISK_STAGE223') for k in row['environment'])
        assert 'Stage223 capability prototype' not in log
    if row['oracle']:
        for marker in ['Block dynamics validation oracle enabled','Block adaptation validation oracle enabled','Block vertical-remap validation oracle enabled']:
            assert marker in log
        assert log.count('remap validation oracle max time')==2
        for k in ['DYNAMICS','ADAPTATION','REMAP']:assert row['environment']['WAVETRISK_VALIDATE_BLOCK_'+k]=='1'
    rows[label]=summary
assert len(r['gates'])==33
for name,g in r['gates'].items():
    assert g['passed'] and g['exact_observed'] and not accept(g['report'],exact=True)
    assert all(v['different']==0 for v in g['report']['fields'].values())
main=subprocess.check_output(['git','rev-parse','main'],text=True).strip();assert main==m['main_commit']
review=dict(stage=226,job='54731',status='passed',archive=str(archive),archive_sha256=digest(archive),inventory_files_verified=len(inventory['files']),source_manifest_sha256=digest(kit/'stage226-sources.json'),builds_verified=6,node=top['host'],placement_counts_verified=list(plans),runs=rows,exact_numerical_gates=33,raw_lifecycle_and_trajectory_rechecked=True,main_commit=main,scope='Restarted J4/J5 resident block pathway on one nk node; no general initialization, arbitrary task-count or multi-node claim; exploratory scaling follows below.',checkpoint_review_limit='Raw checkpoint files excluded from review archive; exact comparison reports re-evaluated and recorded final file hashes cross-checked across runs.')
Path('docs/parallel-blocks-stage226-review.json').write_text(json.dumps(review,indent=2)+'\n')
print('PASS: archive, six build identities, 21 runs, 33 exact reports and raw lifecycle/placement/trajectory evidence')

import statistics as st, csv
from reconcile_evidence_210 import strict_selection, rounding_interval
scaling=[];all_steps=[];baselines={}
for label,row in r['runs'].items():
    if row['kind']!='candidate' or row['oracle']:continue
    assert all(v=='0' for k,v in row['environment'].items() if k.startswith(('WAVETRISK_PROFILE_','WAVETRISK_VALIDATE_')))
    records=analyze(root/label/'run.log')['step_records']
    indices=strict_selection(records)
    assert indices==[2,3,4,5,6,7,10,11,12]
    values=[records[n-1]['seconds'] for n in indices]
    level=row['level'];n=row['ranks'];t=st.median(values)
    baselines.setdefault(level,(n,t,values))
    p0,t0,v0=baselines[level]
    ratios=[b/a for a,b in zip(v0,values)]
    result=dict(level=level,ranks=n,base_ranks=p0,ordinary_indices=indices,ordinary_seconds=values,
        median_seconds=t,stdout_rounding_median_interval=rounding_interval(values,st.median),
        speedup=t0/t,relative_efficiency=t0/t*p0/n,time_change_percent=100*(t/t0-1),
        paired_step_time_ratios=ratios,paired_ratio_range=[min(ratios),max(ratios)],
        median_without_remaps=st.median(records[i-1]['seconds'] for i in indices if not records[i-1]['remap']),
        pre_checkpoint_median=st.median(records[i-1]['seconds'] for i in indices if i<8),
        post_restart_median=st.median(records[i-1]['seconds'] for i in indices if i>9),
        independent_runs=1)
    scaling.append(result)
    for i,step in enumerate(records,1):all_steps.append(dict(run=label,step=i,included=i in indices,**step))
review['scaling']=scaling
review['timing_policy']=dict(included_steps=[2,3,4,5,6,7,10,11,12],excluded_steps=[1,8,9],
    excludes='initial step, checkpoint/save/reload-containing step, first full integration step after that reload',
    remapping='included; removing the two remap steps does not reverse the conclusion',
    clock='rank-zero MPI_Wtime around time_step; stdout cpu label is wall time, not CPU time or measured all-rank maximum',
    precision='stdout ES8.2, three significant digits; the higher-precision run_id_log files were not included in the archive',
    repetition='one independent optimized non-oracle run per task count; the nine steps share run conditions and adaptive history',
    efficiency='E(P;P0)=T(P0)*P0/(T(P)*P), relative to 40 ranks for J4 and 160 for J5, not serial efficiency')
Path('docs/parallel-blocks-stage226-review.json').write_text(json.dumps(review,indent=2)+'\n')
with Path('docs/parallel-blocks-stage226-scaling.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(all_steps[0]));writer.writeheader();writer.writerows(all_steps)
for row in scaling:print(row['level'],row['ranks'],row['median_seconds'],round(row['speedup'],3),round(100*row['relative_efficiency'],1))
