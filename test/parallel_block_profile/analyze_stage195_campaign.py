#!/usr/bin/env python3
"""Audit Stage 195 evidence and reconcile matched phases; no absent-payload recheck."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import statistics as st
from analyze import analyze
from analyze_workflow_195 import collect, summarize, PHASES
from checkpoint_gate import accept
from cluster_194 import workload_gate, stage191_evidence, stage194_evidence
from legacy_metrics_191 import collect as legacy_collect
from rank_metrics_190 import rank_resources
from prepare_cluster_191 import file_digest


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def matched_accounting(runs, mode, select=lambda s: not s['restart']):
    """Use means for additive accounting; never add independently selected maxima."""
    pairs=[]
    for n in range(1,7):
        selected={side:[s for s in runs[f'{side}-{mode}-{n}']['workflow']['steps'] if select(s)]
                  for side in ('block','legacy')}
        b,l=selected.values()
        assert [(s['sequence'],s['restart'],s['remap']) for s in b]==[(s['sequence'],s['restart'],s['remap']) for s in l]
        if not b:continue
        pair={'round':n,'steps':len(b)}
        for stat,wall,phases in [('rank_mean','wall_mean','phase_mean'),('slowest_rank','wall_max','phase_on_slow_rank')]:
            sides={side:{'wall':st.mean(s[wall] for s in rows),
                        'phases':{name:st.mean(s[phases][i] for s in rows) for i,name in enumerate(PHASES)}}
                   for side,rows in selected.items()}
            gap=sides['block']['wall']-sides['legacy']['wall']
            differences={p:sides['block']['phases'][p]-sides['legacy']['phases'][p] for p in PHASES}
            residual=gap-sum(differences.values());assert abs(residual)<2e-6
            pair[stat]={'sides':sides,'gap':gap,'phase_excess':differences,'residual':residual}
        pairs.append(pair)
    if not pairs:return None
    out={'pairs':pairs}
    for stat in ('rank_mean','slowest_rank'):
        out[stat]={'block':st.mean(p[stat]['sides']['block']['wall'] for p in pairs),
                   'legacy':st.mean(p[stat]['sides']['legacy']['wall'] for p in pairs),
                   'gap':st.mean(p[stat]['gap'] for p in pairs),
                   'gap_range':[min(p[stat]['gap'] for p in pairs),max(p[stat]['gap'] for p in pairs)],
                   'phase_excess':{name:st.mean(p[stat]['phase_excess'][name] for p in pairs) for name in PHASES},
                   'max_abs_residual':max(abs(p[stat]['residual']) for p in pairs)}
    return out


def transactions(directory):
    """Ordinary-step rank means by parent phase; nested totals are non-additive."""
    vals=defaultdict(list);counts=defaultdict(list);event_count=0
    intervals=defaultdict(list)
    for file in directory.glob('workflow195-rank-*.txt'):
        flags={}
        for line in file.read_text().splitlines():
            f=line.split()
            if f[0]=='step':flags[int(f[1])]=int(f[2])
            if f[0]=='transaction' and not flags[int(f[1])]:
                key=f'{PHASES[int(f[2])]}:tag{f[3]}'
                vals[key].append(float(f[5]));counts[key].append(int(f[4]))
            if f[0]=='event':
                event_count+=1
                if not flags[int(f[1])] and int(f[2])>0:
                    intervals[(int(file.stem.rsplit('-',1)[1]),int(f[1]),int(f[3]))].append((float(f[5]),float(f[6])))
    # Establish that WT and tendency intervals are disjoint before subtracting both.
    for events in intervals.values():
        events.sort()
        assert all(a[1]<=b[0]+1e-6 for a,b in zip(events,events[1:])), 'Overlapping nested transactions'
    return {'seconds':{k:st.mean(v) for k,v in vals.items()},
            'calls_per_rank_step':{k:st.mean(v) for k,v in counts.items()},'event_count':event_count}


def audit(root, repo, kit):
    root=root.resolve();results=root/'results';r=json.loads((results/'results.json').read_text())
    bpath=root/'builds/wt-stage195/restriction-195-builds.json';b=json.loads(bpath.read_text())
    metadata=json.loads((repo/'docs/parallel-blocks-stage195-kit.json').read_text())
    m=json.loads((kit/'stage195-sources.json').read_text())
    assert sha(Path(str(kit)+'.tar.gz'))==metadata['sha256']
    for name,h in metadata['helpers_sha256'].items():assert sha(kit/name)==h,name
    assert m['candidate_patch_sha256']==metadata['block_patch_sha256']
    assert m['legacy_profile_patch_sha256']==metadata['legacy_probe_patch_sha256']
    assert m['observer_sha256']==metadata['observer_sha256']
    assert m['candidate_sources']['src/workflow_profile_195.f90']==m['legacy_profile_sources']['src/workflow_profile_195.f90']==m['observer_sha256']
    assert b==json.loads(r['build_notes']) and b['status']=='passed'
    assert b['fp_contract']=='off' and b['param']=='param_J5' and b['integrator']=='RK4'
    assert b['baseline_stage']==194 and b['baseline_commit']==metadata['baseline_commit']
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
    assert inventory['stage']==195 and inventory['job']=='54657'
    for name,v in inventory['files'].items():
        path=root/name;assert path.resolve().is_relative_to(root)
        assert path.stat().st_size==v['bytes'] and sha(path)==v['sha256'],name
    assert r['status']=='passed' and len(r['runs'])==58 and len(r['gates'])==146
    assert all(g['passed'] for g in r['gates'].values())
    builds={}
    for kind,v in b['builds'].items():
        path=root/'builds'/Path(v.get('source_root','wt-stage195')).name/f"build-195-{kind.replace('_','-')}.log"
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
    reference=analyze(results/'baseline-warmup/run.log');memory={};tx={};detail={};workflow_files=0
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
            mode=int(env['WAVETRISK_WORKFLOW195']);assert collect(directory,83,mode)==run['workflow'];workflow_files+=83
            if mode>0:tx[label]=transactions(directory)
        else:assert not list(directory.glob('workflow195-rank-*.txt'))
        if run['phase'].startswith('timing-'):
            mode=run['phase'].split('-')[1];assert label==f"{run['side']}-{mode}-{run['round']}"
            assert run['profile']==(mode=='detail')
            if mode!='bare':assert int(env['WAVETRISK_WORKFLOW195'])=={'minimal':0,'coarse':1,'detail':2}[mode]
        if run['profile'] and run['phase']=='timing-detail':
            if run['side']=='block':
                steps=[s for w in a['detail_windows'] for s in w['step_summaries'] if not s['restart']]
                assert len(steps)==11 and all(s['complete'] for s in steps)
                assert max(s['max_self_conservation_error'] for s in steps)<1e-4
                values=[st.mean(s['self_avg'][i] for s in steps) for i in range(67)]
                assert values[56]==0
            else:
                parsed=legacy_collect(directory,83);assert json.loads(json.dumps(parsed))==run['legacy_profile']
                steps=[s for s in parsed['steps'] if not s['restart']];assert len(steps)==11
                values=[st.mean(s['region_mean_self_wall'][i] for s in steps) for i in range(65)]+[0.,0.]
            detail[label]=values
        metrics=rank_resources(directory);assert metrics==run['rank_metrics'] and len(metrics)==83
        assert {v['rank']:v['affinity'] for v in metrics}==placement
        assert all(v['hostname']=='bb13.rhpcs.internal' and v['returncode']==0 for v in metrics)
        assert len(run['memory_before'])==len(run['memory_after'])>0
        for before,after in zip(run['memory_before'],run['memory_after']):
            assert before['host']==after['host']
            assert all(before['vmstat'][k]==after['vmstat'][k] for k in ('pswpin','pswpout'))
        memory[label]={'max_rank_rss_mib':max(v['max_rss_kib'] for v in metrics)/1024,
                       'sum_rank_peak_rss_gib':sum(v['max_rss_kib'] for v in metrics)/1024**2}
    summary=summarize(r['runs'])
    # Python 3.12+ improved float sum: only the cancellation residual differs.
    supplied=json.loads(json.dumps(r['summary']))
    for mode in ('coarse','detail'):
        for got,expected in zip(summary['paired_differences'][mode],supplied['paired_differences'][mode]):
            assert abs(got['residual']-expected['residual'])<1e-12
            expected['residual']=got['residual']
    assert summary==supplied
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
    return {'matched_coarse_excess_ledger':excess_ledger,'status':'passed_archive_audit','job':54657,'inventory_files_verified':len(inventory['files']),
            'workflow_rank_files_reparsed':workflow_files,'results_sha256':sha(results/'results.json'),
            'build_identity_sha256':sha(bpath),'build_logs':builds,'gate_counts':dict(gates),
            'maximum_legacy_field_differences':dict(legacy_max),'summary':summary,'matched_accounting':accounting,
            'transactions':transactions_summary,'transaction_records':tx,'detail_region_self_seconds':detail,
            'memory':memory,'placement_cpus':sorted(c for v in placement.values() for c in v),
            'limitations':['Checkpoints and executables absent: supplied comparison reports and identities audited, payloads not recomputed.',
                           'Matched phase and WT envelopes locate elapsed excess, not its removable fraction or compute/communication cause.',
                           'Tag 2 tendency preparation has different scopes. Never interpret its subtraction as duplicate arithmetic.',
                           'Detailed profiler observer overhead differs by implementation; no uniform correction is valid.',
                           'Local event epochs are not synchronized across ranks. Similar step durations do not exclude waiting or imbalance.',
                           'Mean phase ledger differs from bare median headline; checkpoint steps are separate.',
                           'One 83-rank J5/J7 RK4 fixture on nk; summed peak RSS is not concurrent node RSS.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive_root',type=Path)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--kit',type=Path,required=True)
    p.add_argument('--archive',type=Path,help='Optional original tarball for its digest')
    a=p.parse_args();report=audit(a.archive_root,a.repo,a.kit)
    if a.archive:report['archive_sha256']=sha(a.archive)
    print(json.dumps(report,indent=2))
