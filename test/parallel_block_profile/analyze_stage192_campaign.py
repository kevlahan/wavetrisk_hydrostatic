#!/usr/bin/env python3
"""Audit an extracted Stage 192 review archive without claiming to recheck absent checkpoints."""
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
from cluster_192 import summarize, workload_gate, stage191_evidence, stage192_evidence
from analyze_inverse_192 import compare_profiles
from rank_metrics_190 import rank_resources
from legacy_metrics_191 import collect


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def mean(values): return statistics.mean(values)


TESTED_SOURCE = "30c831149653205210531427f81b184663a3acbf"
TESTED_CHANGED_SOURCES = {
    'src/parallel_block.f90': 'db1d34b7ec1d26c158fc006bc885a8def01ffcda1dbf1d2b3d3d00dba61a5431',
    'src/parallel_block_inverse.f90': 'e3de27f7340af64773f3348317887a8889a150cc2420b11963f72b20a55204e0'}
TESTED_PATCH = 'f23cf5222e285981ffadad1440ec53b23e7d9a74d5eae2436ceed1914f2b6b4c'


def audit(root, repo):
    root=root.resolve();results=root/'results'
    r=json.loads((results/'results.json').read_text())
    build_path=root/'builds/wt-stage192/restriction-192-builds.json'
    b=json.loads(build_path.read_text())
    assert b==json.loads(r['build_notes']) and b['status']=='passed'
    assert b['fp_contract']=='off' and b['param']=='param_J5' and b['integrator']=='RK4'
    assert len(b['physics_objects'])==23
    committed=sources(repo,TESTED_SOURCE)
    assert b['baseline_commit']==TESTED_SOURCE and b['patch_sha256']==TESTED_PATCH
    assert all(committed.get(name)==digest for name,digest in b['baseline_sources'].items())
    committed.update(TESTED_CHANGED_SOURCES)
    assert all(committed.get(name)==digest for name,digest in b['candidate_sources'].items())
    inventory=json.loads((root/'collection.json').read_text())
    assert inventory['stage']==192 and inventory['job']=='54581'
    for name,record in inventory['files'].items():
        path=root/name
        assert path.resolve().is_relative_to(root)
        assert path.stat().st_size==record['bytes'] and sha(path)==record['sha256']
    assert r['status']=='passed' and len(r['runs'])==27 and len(r['gates'])==62
    assert all(g['passed'] for g in r['gates'].values())
    summary=summarize(r['runs']);assert summary==r['summary']
    logs={}
    for kind,record in b['builds'].items():
        directory=root/'builds'/Path(record.get('source_root','wt-stage192')).name
        name='legacy-profile' if kind=='legacy_profile' else kind
        path=directory/f'build-192-{name}.log'
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
        assert a['non_checkpoint_median']==run['ordinary_step_median']
        text=(results/label/'run.log').read_text()
        import re
        assert list(map(int,re.findall(r'Restarting from checkpoint\s+(\d+)',text)))==[3,4]
        assert list(map(int,re.findall(r'Saving checkpoint\s+(\d+)',text)))==[4]
        assert 'Remapping vertical coordinates' in text.split('Restarting from checkpoint',2)[2]
        env=run['environment']
        assert env['WAVETRISK_STAGE191_INSTALL']=='1' and env['WAVETRISK_STAGE191_MASS_OVERLAP']=='0'
        if run['kind'] in ('candidate','checked'):
            sparse='0' if label=='candidate-control' else '1'
            assert stage192_evidence(results/label,sparse)==run['stage192_evidence']
            assert stage191_evidence(results/label,a,'1','0',run['profile'])==run['stage191_evidence']
        for name in ('DYNAMICS','ADAPTATION','REMAP'):
            assert env['WAVETRISK_VALIDATE_BLOCK_'+name]==('1' if run['oracles'] else '0')
        metrics=run['rank_metrics'];assert len(metrics)==83
        assert rank_resources(results/label)==metrics
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
        if kind=='candidate':assert stage191_evidence(results/'candidate-profile',a,'1','0',True)==r['runs']['candidate-profile']['stage191_evidence']
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
    # Shared region IDs allow before/after block attribution; keep it separate
    # from legacy's differently nested component timers.
    components={}
    for kind in ('baseline','candidate'):
        parsed=analyze(results/(kind+'-profile')/'run.log')
        components[kind]={name:sum(v['self_wall_avg'] for w in parsed['detail_windows'] for v in w['regions'] if v['id']==i)
            for i,name in [(15,'plan'),(16,'pack'),(17,'local_aliases'),(18,'MPI_post_self_copy'),(19,'MPI_wait'),(20,'install')]}
    matched=compare_profiles(results,83)
    assert matched==r['inverse_comparison']
    assert matched==json.loads((results/'inverse-comparison.json').read_text())
    return {'job':54581,'status':'passed_review','inventory_files_verified':len(inventory['files']),
        'source_hashes_match_frozen_candidate':True,'baseline_commit':TESTED_SOURCE,'candidate_patch_sha256':TESTED_PATCH,
        'matched_inverse_comparison':matched,'inverse_component_rank_mean_self_seconds':components,
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
        'Only one profile per implementation. Profile component differences are diagnostic, not repeated performance measurements.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('extracted_archive',type=Path);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
    print(json.dumps(audit(a.extracted_archive,a.repo),indent=2))
