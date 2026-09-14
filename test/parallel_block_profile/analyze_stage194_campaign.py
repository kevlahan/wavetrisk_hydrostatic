#!/usr/bin/env python3
"""Audit an extracted Stage 194 review archive without claiming to recheck absent checkpoints."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
import re
from analyze import analyze
from analyze_stage190_gap import GROUPS
from prepare_cluster_191 import file_digest
from checkpoint_gate import accept
from cluster_194 import summarize, workload_gate, stage191_evidence, stage194_evidence
from analyze_inner_inverse_194 import compare_profiles
from rank_metrics_190 import rank_resources
from legacy_metrics_191 import collect


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def mean(values): return statistics.mean(values)


MANIFEST_SHA256 = '751abec0eae507384d19dd84e20181f3908df03010587051a28dbc39daa612d2'


def audit(root, repo, kit):
    root=root.resolve();results=root/'results'
    r=json.loads((results/'results.json').read_text())
    build_path=root/'builds/wt-stage194/restriction-194-builds.json'
    b=json.loads(build_path.read_text())
    assert b==json.loads(r['build_notes']) and b['status']=='passed'
    assert b['fp_contract']=='off' and b['param']=='param_J5' and b['integrator']=='RK4'
    assert len(b['physics_objects'])==23
    assert b['build_host']['node']=='bb13.rhpcs.internal'
    manifest_path=kit/'stage194-sources.json'
    assert sha(manifest_path)==MANIFEST_SHA256
    m=json.loads(manifest_path.read_text())
    for key in ('baseline_commit','baseline_stage','baseline_patch_sha256','baseline_sources','candidate_sources','patch_sha256','makefile_sha256'):
        assert b[key]==m['candidate_patch_sha256' if key=='patch_sha256' else key],key
    assert b['baseline_stage']==193
    legacy=b['builds']['legacy'];lp=b['builds']['legacy_profile']
    assert legacy['sources']==m['legacy_io_sources']
    assert legacy['original_sources']==m['legacy_original_sources']
    assert legacy['io_patch_sha256']==m['legacy_io_patch_sha256']
    assert lp['sources']==m['legacy_profile_sources']
    for n,h in m['candidate_sources'].items():
        assert file_digest(repo/n)==h,n
    inventory=json.loads((root/'collection.json').read_text())
    assert inventory['stage']==194 and inventory['job']=='54653'
    for name,record in inventory['files'].items():
        path=root/name
        assert path.resolve().is_relative_to(root)
        assert path.stat().st_size==record['bytes'] and sha(path)==record['sha256']
    assert r['status']=='passed' and len(r['runs'])==28 and len(r['gates'])==64
    assert all(g['passed'] for g in r['gates'].values())
    summary=summarize(r['runs']);assert summary==r['summary']
    logs={}
    for kind,record in b['builds'].items():
        directory=root/'builds'/Path(record.get('source_root','wt-stage194')).name
        name='legacy-profile' if kind=='legacy_profile' else kind
        path=directory/f'build-194-{name}.log'
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
    assert len(placement)==83 and all(len(v)==1 for v in placement.values()) and len({c for v in placement.values() for c in v})==83
    memory={};reference=analyze(results/'baseline-warmup/run.log')
    for label,run in r['runs'].items():
        assert run['returncode']==0 and run['binary_sha256']==b['builds'][run['kind']]['sha256']
        a=analyze(results/label/'run.log');assert a['completed']
        assert workload_gate(reference,a)['passed']
        assert len(a['step_records'])==run['step_count']==12
        assert run['reloads']==[3,4] and run['post_restart_remaps']['4']
        assert a['non_checkpoint_median']==run['ordinary_step_median']
        text=(results/label/'run.log').read_text()
        assert list(map(int,re.findall(r'Restarting from checkpoint\s+(\d+)',text)))==[3,4]
        assert list(map(int,re.findall(r'Saving checkpoint\s+(\d+)',text)))==[4]
        assert 'Remapping vertical coordinates' in text.split('Restarting from checkpoint',2)[2]
        env=run['environment']
        assert env['WAVETRISK_STAGE191_INSTALL']=='1' and env['WAVETRISK_STAGE191_MASS_OVERLAP']=='0'
        assert env['WAVETRISK_STAGE192_SPARSE_OUTER']=='0'
        assert env['WAVETRISK_STAGE193_SCALAR_COLUMNS']=='1'
        if run['kind'] in ('baseline','candidate','checked'):
            assert re.findall(r'Stage 193 requested scalar inverse columns = ([01])',text)==['1']
            assert re.findall(r'Stage 192 requested sparse outer install = ([01])',text)==['0']
            evidence=stage191_evidence(results/label,a,'1','0',run['profile'])
            if 'stage191_evidence' in run: assert evidence==run['stage191_evidence']
        if run['kind'] in ('candidate','checked'):
            setting='0' if label.startswith('candidate-control') else '1'
            assert env['WAVETRISK_STAGE194_INNER_GEOMETRY']==setting
            assert stage194_evidence(results/label,setting)==run['stage194_evidence']
        oracle_messages=('Block dynamics validation oracle enabled','Block adaptation validation oracle enabled',
                         'Block vertical-remap validation oracle enabled')
        if run['kind'] in ('candidate','checked'):
            assert all((message in text)==run['oracles'] for message in oracle_messages)
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
    for kind in ('baseline','candidate-control','candidate'):
        a=analyze(results/(kind+'-profile')/'run.log')
        steps=[s for w in a['detail_windows'] for s in w['step_summaries'] if not s['restart']]
        assert len(steps)==11 and all(s['complete'] for s in steps)
        error=max(s['max_self_conservation_error'] for s in steps);assert error<1e-4
        values=[mean(s['self_avg'][i] for s in steps) for i in range(len(steps[0]['self_avg']))]
        values += [0.]*(67-len(values))
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
        'mean_self_wall_by_region':[mean(s['region_mean_self_wall'][i] for s in steps) for i in range(65)]+[0.,0.]}
    groups_ids={k:list(v) for k,v in GROUPS.items()}
    groups_ids['Native inverse'] += [66,67]
    assert sorted(i for ids in groups_ids.values() for i in ids)==list(range(1,68))
    groups={name:{kind:sum(v['mean_self_wall_by_region'][i-1] for i in ids) for kind,v in profiles.items()} for name,ids in groups_ids.items()}
    for kind,profile in profiles.items():
        assert abs(sum(g[kind] for g in groups.values())-profile['mean_rank_mean_wall'])<1e-4
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
    pairs=[]
    for i in range(1,7):
        base=r['runs'][f'baseline-{i}'];candidate=r['runs'][f'candidate-{i}']
        pairs.append({'round':i,'elapsed_saved_seconds':base['elapsed_seconds']-candidate['elapsed_seconds'],
                      'ordinary_step_saved_seconds':base['ordinary_step_median']-candidate['ordinary_step_median']})
    return {'paired_timings':pairs,'candidate_control':{k:r['runs']['candidate-control'][k] for k in
            ('elapsed_seconds','ordinary_step_median','stage194_evidence')},'job':54653,'status':'passed_review','inventory_files_verified':len(inventory['files']),
        'source_hashes_match_frozen_candidate':True,'baseline_commit':b['baseline_commit'],'baseline_stage':193,'candidate_patch_sha256':b['patch_sha256'],
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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('extracted_archive',type=Path);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--kit',type=Path,required=True);a=p.parse_args()
    print(json.dumps(audit(a.extracted_archive,a.repo,a.kit),indent=2))
