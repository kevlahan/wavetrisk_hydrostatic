#!/usr/bin/env python3
"""Audit Stage 196 evidence and reconcile matched phases; no absent-payload recheck."""
import argparse
from collections import Counter
import hashlib
import math
import json
from pathlib import Path
import re
import statistics as st
from analyze import analyze
from analyze_workflow_195 import collect
from checkpoint_gate import accept
from cluster_194 import workload_gate, stage191_evidence, stage194_evidence
from analyze_stage195_campaign import matched_accounting, transactions
from analyze_operations_196 import collect as operation_collect, summarize as operation_summary, GROUPS, BOUNDARY_GROUPS
from cluster_196 import operation_mode
from rank_metrics_190 import rank_resources
from prepare_cluster_191 import file_digest


def close(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a:close(a[k],b[k],path+'/'+str(k))
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):close(x,y,path+'/'+str(i))
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-12),(path,a,b)
    else:assert a==b,(path,a,b)


def operation_accounting(runs,summary):
    out={}
    for mode,pairs in summary['operation_excess'].items():
        contexts={}
        for c in ('WT','outside_WT'):
            groups={g:{side:st.mean(p['contexts'][c]['sides'][side][g] for p in pairs) for side in ('block','legacy')} for g in GROUPS}
            for v in groups.values():v['excess']=v['block']-v['legacy']
            contexts[c]={'groups':groups,'excess':st.mean(p['contexts'][c]['total_excess'] for p in pairs),
                         'boundary_transfer_excess':st.mean(p['contexts'][c]['boundary_and_transfer_excess'] for p in pairs)}
            for side in ('block','legacy'):
                rows=[step['contexts'][c] for n in range(1,7) for step in runs[f'{side}-{mode}-{n}']['operations']['steps'] if not step['restart']]
                contexts[c][side]={'self_seconds':[st.mean(s['self_seconds'][i] for s in rows) for i in range(78)],
                    'calls_per_rank_step':[st.mean(s['calls_global'][i]/83 for s in rows) for i in range(78)],
                    'counter_per_rank_step':[st.mean(s['counters_global'][i]/83 for s in rows) for i in range(27)]}
        out[mode]=contexts
    return out


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def audit(root, repo, kit):
    root=root.resolve();results=root/'results';r=json.loads((results/'results.json').read_text())
    bpath=root/'builds/wt-stage196/restriction-196-builds.json';b=json.loads(bpath.read_text())
    metadata=json.loads((repo/'docs/parallel-blocks-stage196-kit.json').read_text())
    m=json.loads((kit/'stage196-sources.json').read_text())
    assert sha(Path(str(kit)+'.tar.gz'))==metadata['sha256']
    for name,h in metadata['helpers_sha256'].items():assert sha(kit/name)==h,name
    assert m['candidate_patch_sha256']==metadata['candidate_patch_sha256']
    assert m['legacy_profile_patch_sha256']==metadata['legacy_profile_patch_sha256']
    assert m['observer_sha256']==metadata['observer_sha256']
    assert m['candidate_sources']['src/workflow_profile_195.f90']==m['legacy_profile_sources']['src/workflow_profile_195.f90']==m['observer_sha256']
    assert b==json.loads(r['build_notes']) and b['status']=='passed'
    assert b['fp_contract']=='off' and b['param']=='param_J5' and b['integrator']=='RK4'
    assert b['baseline_stage']==195 and b['baseline_commit']==metadata['baseline_commit']
    assert len(b['physics_objects'])==23 and b['build_host']['node']=='bb13.rhpcs.internal'
    for key in ('baseline_commit','baseline_stage','baseline_sources','candidate_sources','patch_sha256','makefile_sha256'):
        assert b[key]==m['candidate_patch_sha256' if key=='patch_sha256' else key],key
    assert b['builds']['legacy']['sources']==m['legacy_io_sources']
    assert b['builds']['legacy']['original_sources']==m['legacy_original_sources']
    assert b['builds']['legacy']['io_patch_sha256']==m['legacy_io_patch_sha256']
    assert b['builds']['legacy']['revision']==m['legacy_revision']
    assert b['builds']['legacy_profile']['sources']==m['legacy_profile_sources']
    for name,h in m['baseline_sources'].items():assert file_digest(repo/name)==h,name
    inventory=json.loads((root/'collection.json').read_text())
    assert inventory['stage']==196 and inventory['job']=='54659'
    for name,v in inventory['files'].items():
        path=root/name;assert path.resolve().is_relative_to(root)
        assert path.stat().st_size==v['bytes'] and sha(path)==v['sha256'],name
    assert r['status']=='passed' and len(r['runs'])==58 and len(r['gates'])==146
    assert all(g['passed'] for g in r['gates'].values())
    builds={}
    for kind,v in b['builds'].items():
        path=root/'builds'/Path(v.get('source_root','wt-stage196')).name/f"build-196-{kind.replace('_','-')}.log"
        assert sha(path)==v['log_sha256']
        lines=[l for l in path.read_text().splitlines() if ' -c src/' in l and 'mpif90 ' in l]
        assert lines and all('-ffp-contract=off' in l for l in lines)
        if kind=='checked':assert all(' -O0 ' in l and '-fcheck=all' in l for l in lines)
        else:assert all(' -O2 ' in l and '-march=native' in l for l in lines)
        builds[kind]={'sha256':sha(path),'compiler_commands_verified':len(lines)}
    gates=Counter();legacy_max=Counter()
    for name,g in r['gates'].items():
        if name.startswith('workload-'):gates['workload']+=1;continue
        exact=not name.startswith('legacy-screen-');gates['exact' if exact else 'legacy_screen']+=1
        for report in g['checkpoints'].values():
            assert not accept(report,exact=exact)
            if not exact:
                for field,v in report['fields'].items():legacy_max[field.split('/')[1]]=max(legacy_max[field.split('/')[1]],v['max_abs'])
    placement={v['rank']:v['allowed_cpus'] for v in r['placement']['ranks']}
    assert len(placement)==83 and all(len(v)==1 for v in placement.values())
    assert len({c for v in placement.values() for c in v})==83
    reference=analyze(results/'baseline-warmup/run.log');memory={};tx={};workflow_files=0;operation_files=0
    for label,run in r['runs'].items():
        directory=results/label;a=analyze(directory/'run.log');text=(directory/'run.log').read_text()
        assert run['returncode']==0 and run['binary_sha256']==b['builds'][run['kind']]['sha256']
        assert a['completed'] and workload_gate(reference,a)['passed']
        assert len(a['step_records'])==run['step_count']==12
        assert run['ordinary_step_median']==a['non_checkpoint_median']
        assert run['reloads']==[3,4] and run['post_restart_remaps']['4']
        assert list(map(int,re.findall(r'Restarting from checkpoint\s+(\d+)',text)))==[3,4]
        assert list(map(int,re.findall(r'Saving checkpoint\s+(\d+)',text)))==[4]
        assert 'Remapping vertical coordinates' in text.split('Restarting from checkpoint',2)[2]
        env=run['environment']
        for key,val in [('STAGE191_INSTALL','1'),('STAGE191_MASS_OVERLAP','0'),('STAGE192_SPARSE_OUTER','0'),('STAGE193_SCALAR_COLUMNS','1'),('STAGE194_INNER_GEOMETRY','1')]:assert env['WAVETRISK_'+key]==val
        for name in ('DYNAMICS','ADAPTATION','REMAP'):assert env['WAVETRISK_VALIDATE_BLOCK_'+name]==str(int(run['oracles']))
        if run['side']=='block':
            assert re.findall(r'Stage 193 requested scalar inverse columns = ([01])',text)==['1']
            assert re.findall(r'Stage 192 requested sparse outer install = ([01])',text)==['0']
            if 'stage191_evidence' in run:assert stage191_evidence(directory,a,'1','0',run['profile'])==run['stage191_evidence']
            if 'stage194_evidence' in run:assert stage194_evidence(directory,'1')==run['stage194_evidence']
            for message in ('Block dynamics validation oracle enabled','Block adaptation validation oracle enabled','Block vertical-remap validation oracle enabled'):assert (message in text)==run['oracles']
        if 'workflow' in run:
            mode=int(env['WAVETRISK_WORKFLOW195']);close(collect(directory,83,mode),run['workflow'],label+'/workflow');workflow_files+=83
            if mode>0:tx[label]=transactions(directory)
        else:assert not list(directory.glob('workflow195-rank-*.txt'))
        assert not run['profile']
        for key in ('WAVETRISK_PROFILE_BLOCK_DETAIL','WAVETRISK_PROFILE_PARALLEL_BLOCKS'):
            assert env[key]=='0'
        if 'workflow' in run:
            assert int(env['WAVETRISK_WORKFLOW195'])==1
            mode=operation_mode(label,run['phase'])
            assert int(env['WAVETRISK_OPERATIONS196'])==mode
            parsed=operation_collect(directory,83,mode,run['side'])
            close(parsed,run['operations'],label+'/operations')
            operation_files+=83
        else:assert not list(directory.glob('operations196-rank-*.txt'))
        if run['phase'].startswith('timing-'):
            mode=run['phase'].split('-')[1]
            assert label==f"{run['side']}-{mode}-{run['round']}"
        metrics=rank_resources(directory);assert metrics==run['rank_metrics'] and len(metrics)==83
        assert {v['rank']:v['affinity'] for v in metrics}==placement
        assert all(v['hostname']=='bb13.rhpcs.internal' and v['returncode']==0 for v in metrics)
        assert len(run['memory_before'])==len(run['memory_after'])>0
        for before,after in zip(run['memory_before'],run['memory_after']):
            assert before['host']==after['host']
            assert all(before['vmstat'][k]==after['vmstat'][k] for k in ('pswpin','pswpout'))
        memory[label]={'max_rank_rss_mib':max(v['max_rss_kib'] for v in metrics)/1024,
                       'sum_rank_peak_rss_gib':sum(v['max_rss_kib'] for v in metrics)/1024**2}
    summary=operation_summary(r['runs'])
    close(summary,r['summary'],'summary')
    accounting={mode:{name:matched_accounting(r['runs'],mode,select) for name,select in
                [('ordinary',lambda s:not s['restart']),('ordinary_without_remap',lambda s:not s['restart'] and not s['remap']),
                 ('ordinary_with_remap',lambda s:not s['restart'] and bool(s['remap'])),('checkpoint',lambda s:bool(s['restart']))]}
                for mode in ('coarse','detail')}
    transactions_summary={mode:{side:{key:st.mean(tx[f'{side}-{mode}-{n}']['seconds'][key] for n in range(1,7))
                         for key in tx[f'{side}-{mode}-1']['seconds']} for side in ('block','legacy')} for mode in ('coarse','detail')}
    coarse=accounting['coarse']['ordinary']['rank_mean']
    t=transactions_summary['coarse']
    wt=t['block']['advance:tag1']-t['legacy']['advance:tag1']
    excess_ledger={'dynamics_WT_after_step':wt,
                  'dynamics_outside_WT_after_step':coarse['phase_excess']['advance']-wt,
                  'adaptation_reconstruction':coarse['phase_excess']['adapt_reconstruct'],
                  'other_phases':coarse['gap']-coarse['phase_excess']['advance']-coarse['phase_excess']['adapt_reconstruct']}
    assert abs(sum(excess_ledger.values())-coarse['gap'])<1e-12
    counter_limits={'unavailable_ids':[4,5,12,13,14,15],
                    'reason':'Call-site detail_enabled guards suppress these counters with the old profiler off. Zero is unavailable, not zero work.',
                    'source_locations':['src/comm_mpi.f90:1037','src/comm_mpi.f90:1066','src/comm_mpi.f90:1094','src/comm_mpi.f90:1510','src/parallel_block_mpi.f90:20980']}
    for run in r['runs'].values():
        if 'operations' not in run:continue
        for step in run['operations']['steps']:
            for context in ('WT','outside_WT'):
                assert all(step['contexts'][context]['counters_global'][i-1]==0 for i in counter_limits['unavailable_ids'])
    return {'counter_limitations':counter_limits,'operation_accounting':operation_accounting(r['runs'],summary),'matched_coarse_excess_ledger':excess_ledger,'status':'passed_archive_audit','job':54659,'inventory_files_verified':len(inventory['files']),
            'operation_rank_files_reparsed':operation_files,'workflow_rank_files_reparsed':workflow_files,'results_sha256':sha(results/'results.json'),
            'build_identity_sha256':sha(bpath),'build_logs':builds,'gate_counts':dict(gates),
            'maximum_legacy_field_differences':dict(legacy_max),'summary':summary,'matched_accounting':accounting,
            'transactions':transactions_summary,'transaction_records':tx,
            'memory':memory,'placement_cpus':sorted(c for v in placement.values() for c in v),
            'limitations':['Checkpoints and executables absent: supplied comparison reports and identities audited, payloads not recomputed.',
                           'Matched operation ledgers locate elapsed excess; MPI wait may reflect late computation/progress, not network service time.',
                           'Tag 2 tendency preparation has different scopes. Never interpret its subtraction as duplicate arithmetic.',
                           'Operation observer effects calibrated against minimal; minimal itself has a legacy cost relative to bare. No uniform correction applied.',
                           'Local event epochs are not synchronized across ranks. Similar step durations do not exclude waiting or imbalance.',
                           'Mean phase ledger differs from bare median headline; checkpoint steps are separate.',
                           'One 83-rank J5/J7 RK4 fixture on nk; summed peak RSS is not concurrent node RSS.',
                           'Counters 4,5,12,13,14,15 suppressed by disabled old-profiler call-site guards; zero is unavailable.',
                           'Kernel groups include different loop/setup scopes across implementations; their differences are not per-operation arithmetic slowdowns.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive_root',type=Path)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--kit',type=Path,required=True)
    p.add_argument('--archive',type=Path,help='Optional original tarball for its digest')
    a=p.parse_args();report=audit(a.archive_root,a.repo,a.kit)
    if a.archive:report['archive_sha256']=sha(a.archive)
    print(json.dumps(report,indent=2))
