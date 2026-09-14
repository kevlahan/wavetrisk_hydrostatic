#!/usr/bin/env python3
"""Audit an extracted Stage 191 review archive without claiming to recheck absent checkpoints."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
from analyze import analyze
from analyze_stage190_gap import GROUPS
from prepare_cluster_191 import sources
from checkpoint_gate import accept
from cluster_191 import summarize, workload_gate, stage191_evidence
from legacy_metrics_191 import collect


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def mean(values): return statistics.mean(values)


TESTED_SOURCE = "37ac536405afbb8b2f8d59037db126950af574cd"


def audit(root, repo):
    root=root.resolve();results=root/'results'
    r=json.loads((results/'results.json').read_text())
    build_path=root/'builds/wt-stage191/restriction-191-builds.json'
    b=json.loads(build_path.read_text())
    assert b==json.loads(r['build_notes']) and b['status']=='passed'
    assert b['fp_contract']=='off' and b['param']=='param_J5' and b['integrator']=='RK4'
    assert len(b['physics_objects'])==23
    committed=sources(repo,TESTED_SOURCE)
    assert all(committed.get(name)==digest for name,digest in b['candidate_sources'].items())
    assert r['status']=='passed' and len(r['runs'])==33 and len(r['gates'])==74
    assert all(g['passed'] for g in r['gates'].values())
    summary=summarize(r['runs']);assert summary==r['summary']
    logs={}
    for kind,record in b['builds'].items():
        directory=root/'builds'/Path(record.get('source_root','wt-stage191')).name
        name='legacy-profile' if kind=='legacy_profile' else kind
        path=directory/f'build-191-{name}.log'
        assert sha(path)==record['log_sha256']
        compiler_lines=[line for line in path.read_text().splitlines() if ' -c src/' in line and 'mpif90 ' in line]
        assert compiler_lines and all('-ffp-contract=off' in line for line in compiler_lines)
        if kind=='checked':assert all(' -O0 ' in line and '-fcheck=all' in line for line in compiler_lines)
        else:assert all(' -O2 ' in line and '-march=native' in line for line in compiler_lines)
        logs[kind]={'sha256':sha(path),'compiler_commands_verified':len(compiler_lines)}
    counts=Counter();legacy_max=Counter()
    for label,gate in r['gates'].items():
        if label.startswith('workload-'): counts['workload']+=1;continue
        exact=not label.startswith('legacy-screen-');counts['exact' if exact else 'legacy_screen']+=1
        for report in gate['checkpoints'].values():
            assert not accept(report,exact=exact)
            if not exact:
                for field,v in report['fields'].items():
                    key=field.split('/')[1];legacy_max[key]=max(legacy_max[key],v['max_abs'])
    placement={v['rank']:v['allowed_cpus'] for v in r['placement']['ranks']}
    assert len(placement)==83 and len({c for v in placement.values() for c in v})==83
    memory={};reference=analyze(results/'baseline-warmup/run.log')
    for label,run in r['runs'].items():
        assert run['returncode']==0 and run['binary_sha256']==b['builds'][run['kind']]['sha256']
        a=analyze(results/label/'run.log');assert a['completed']
        assert workload_gate(reference,a)['passed']
        assert len(a['step_records'])==run['step_count']==12
        assert run['reloads']==[3,4] and run['post_restart_remaps']['4']
        metrics=run['rank_metrics'];assert len(metrics)==83
        assert {v['rank']:v['affinity'] for v in metrics}==placement
        assert all(v['hostname']=='bb13.rhpcs.internal' and v['returncode']==0 for v in metrics)
        delta=[]
        for before,after in zip(run['memory_before'],run['memory_after']):
            assert before['host']==after['host']
            delta.append({k:after['vmstat'][k]-before['vmstat'][k] for k in ('pswpin','pswpout')})
        assert delta and all(v=={'pswpin':0,'pswpout':0} for v in delta)
        memory[label]={'max_rank_rss_mib':max(v['max_rss_kib'] for v in metrics)/1024,
                       'sum_rank_peak_rss_gib':sum(v['max_rss_kib'] for v in metrics)/1024**2,
                       'major_faults_sum':sum(v['major_faults'] for v in metrics),'swap_deltas':delta}
    profiles={}
    for kind in ('baseline','candidate'):
        a=analyze(results/(kind+'-profile')/'run.log')
        steps=[s for w in a['detail_windows'] for s in w['step_summaries'] if not s['restart']]
        assert len(steps)==11 and all(s['complete'] for s in steps)
        error=max(s['max_self_conservation_error'] for s in steps);assert error<1e-4
        values=[mean(s['self_avg'][i] for s in steps) for i in range(65)]
        counters=Counter()
        for w in a['detail_windows']:counters.update({k.strip():v for k,v in w['counters'].items()})
        profiles[kind]={'ordinary_steps':[s['sequence'] for s in steps],
            'mean_rank_mean_wall':mean(s['wall_avg'] for s in steps),'mean_self_wall_by_region':values,
            'max_conservation_error':error,'whole_window_counters':dict(counters)}
        if kind=='candidate':assert stage191_evidence(results/'candidate-profile',a,'1','1',True)==r['runs']['candidate-profile']['stage191_evidence']
    legacy=collect(results/'legacy-profile',83)
    assert json.loads(json.dumps(legacy))==r['runs']['legacy-profile']['legacy_profile']
    steps=[s for s in legacy['steps'] if not s['restart']];assert len(steps)==11
    profiles['legacy']={'ordinary_steps':[s['sequence'] for s in steps],
        'mean_rank_mean_wall':mean(s['wall_mean'] for s in steps),
        'mean_self_wall_by_region':[mean(s['region_mean_self_wall'][i] for s in steps) for i in range(65)]}
    groups={name:{kind:sum(v['mean_self_wall_by_region'][i-1] for i in ids) for kind,v in profiles.items()} for name,ids in GROUPS.items()}
    rss={}
    for kind in ('baseline','candidate','legacy'):
        selected=[memory[k] for k,v in r['runs'].items() if v['kind']==kind and v['phase']=='timing']
        rss[kind]={'max_rank_rss_mib_range':[min(v['max_rank_rss_mib'] for v in selected),max(v['max_rank_rss_mib'] for v in selected)],
                   'median_sum_rank_peak_rss_gib':statistics.median(v['sum_rank_peak_rss_gib'] for v in selected)}
    return {'job':54561,'status':'passed_review','source_hashes_match_committed_candidate':True,'source_commit':TESTED_SOURCE,
        'results_sha256':sha(results/'results.json'),'build_identity_sha256':sha(build_path),
        'build_logs':logs,'gate_counts':dict(counts),'maximum_legacy_field_differences':dict(legacy_max),
        'summary':summary,'ordinary_step_reduction_percent':100*(1-summary['candidate']['median_ordinary_step_seconds']/summary['baseline']['median_ordinary_step_seconds']),
        'fraction_of_excess_over_legacy_removed_percent':100*(summary['baseline']['median_ordinary_step_seconds']-summary['candidate']['median_ordinary_step_seconds'])/(summary['baseline']['median_ordinary_step_seconds']-summary['legacy']['median_ordinary_step_seconds']),
        'memory':memory,'timing_memory_summary':rss,'profiles':profiles,'profile_groups':groups,
        'limitations':['Checkpoint payloads and executables are absent: validate supplied reports/hashes, not re-execution.',
        'One short 83-rank nk J5/J7 RK4 campaign; no general multi-node or longer-run claim.',
        'Legacy phase boundaries differ, notably boundary work split out of inverse/scalar phases. Group differences are not a causal excess decomposition.',
        'Profile means exclude checkpoint steps but are instrumented runs; timing medians exclude all profiles.',
        'Summed rank peaks are not simultaneous node memory; resource collection excludes waiting Python parents.',
        'Ablations follow the rotating timing triples; overlap has no demonstrated gain here, not proof it can never help.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('extracted_archive',type=Path);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
    print(json.dumps(audit(a.extracted_archive,a.repo),indent=2))
