#!/usr/bin/env python3
"""Account for paired legacy/block profile time without treating labels as causes."""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics
from analyze import analyze
from analyze_stage190_gap import GROUPS
from analyze_inverse_192 import compare_profiles
from legacy_metrics_191 import collect

GROUP_IDS={k:list(v) for k,v in GROUPS.items()}
GROUP_IDS['Native inverse'] += [66,67]
assert sorted(i for ids in GROUP_IDS.values() for i in ids)==list(range(1,68))


def check_records(records,ranks):
    if not records: raise ValueError('No step records')
    for seq,rows in records.items():
        if set(rows)!=set(range(ranks)): raise ValueError('Incomplete rank coverage')
        if len({(r['restart'],r['remap']) for r in rows.values()})!=1:
            raise ValueError('Rank event flags disagree')
        for r in rows.values():
            if len(r['self'])!=67: raise ValueError('Wrong region count')
            if any(not math.isfinite(x) or x<0 for x in [r['wall'],*r['self']]):
                raise ValueError('Invalid time')
            if abs(sum(r['self'])-r['wall'])>1e-4: raise ValueError('Step self times do not conserve wall')


def block_records(path,ranks):
    parsed=analyze(path)
    if not parsed['completed']: raise ValueError('Incomplete block log')
    records=defaultdict(dict)
    for w in parsed['detail_windows']:
        for s in w['steps'].values():
            # Parser keys carry rank and sequence; records also expose both.
            seq=s['sequence'];rank=s['rank']
            if rank in records[seq]: raise ValueError('Duplicate block rank/step')
            values=s['self-wall']
            if len(values) not in (65,67): raise ValueError('Unsupported block schema')
            records[seq][rank]={'wall':s['wall'],'self':values+[0.]*(67-len(values)),
                                'restart':bool(s['restarts']),'remap':bool(s['remaps'])}
    check_records(records,ranks)
    return records


def legacy_records(path,ranks):
    collect(path,ranks)  # Check names, records and standard schema first.
    records=defaultdict(dict)
    for file in path.glob('legacy-detail-rank-*.txt'):
        rank=int(file.stem.rsplit('-',1)[1]);ordinal={}
        for line in file.read_text().splitlines():
            f=line.split()
            if not f: continue
            if f[0]=='step':
                seq=int(f[1]);ordinal[len(ordinal)+1]=seq
                records[seq][rank]={'wall':float(f[4]),'self':[0.]*67,
                                    'restart':bool(int(f[2])),'remap':bool(int(f[3]))}
            elif f[0]=='self':
                records[ordinal[int(f[1])]][rank]['self'][int(f[2])-1]=float(f[3])
    check_records(records,ranks)
    return records


def ledger(block,legacy):
    """Signed exclusive accounting, using each implementation's actual slow rank."""
    if set(block)!=set(legacy): raise ValueError('Unmatched step sequence')
    rows=[]
    for seq in sorted(block):
        b,l=block[seq],legacy[seq]
        if set(b)!=set(l): raise ValueError('Unmatched ranks')
        if (b[0]['restart'],b[0]['remap'])!=(l[0]['restart'],l[0]['remap']):
            raise ValueError('Unmatched restart/remap workload')
        if b[0]['restart']:continue
        entry={'sequence':seq,'remap':b[0]['remap']}
        for mode in ('rank_mean','slowest_rank'):
            values={};walls={};chosen={}
            for kind,rec in [('block',b),('legacy',l)]:
                slow=max(rec,key=lambda r:rec[r]['wall']);chosen[kind]=slow
                if mode=='rank_mean':
                    walls[kind]=statistics.mean(r['wall'] for r in rec.values())
                    values[kind]=[statistics.mean(r['self'][i] for r in rec.values()) for i in range(67)]
                else:
                    walls[kind]=rec[slow]['wall'];values[kind]=rec[slow]['self']
            groups={name:{k:sum(values[k][i-1] for i in ids) for k in values} for name,ids in GROUP_IDS.items()}
            for g in groups.values():g['difference_seconds']=g['block']-g['legacy']
            gap=walls['block']-walls['legacy']
            residual=gap-sum(g['difference_seconds'] for g in groups.values())
            if abs(residual)>2e-4: raise ValueError('Signed group differences fail conservation')
            entry[mode]={'wall':walls,'gap_seconds':gap,'groups':groups,'rounding_residual_seconds':residual}
            if mode=='slowest_rank':entry[mode]['selected_ranks']=chosen
        rows.append(entry)
    if not rows: raise ValueError('No ordinary paired steps')
    summaries={}
    for label,selected in [('all_ordinary',rows),('without_remap',[r for r in rows if not r['remap']]),
                           ('with_remap',[r for r in rows if r['remap']])]:
        if not selected:continue
        summaries[label]={'steps':len(selected)}
        for mode in ('rank_mean','slowest_rank'):
            summaries[label][mode]={
                'block_seconds':statistics.mean(r[mode]['wall']['block'] for r in selected),
                'legacy_seconds':statistics.mean(r[mode]['wall']['legacy'] for r in selected),
                'gap_seconds':statistics.mean(r[mode]['gap_seconds'] for r in selected),
                'group_differences_seconds':{n:statistics.mean(r[mode]['groups'][n]['difference_seconds'] for r in selected) for n in GROUP_IDS},
                'rounding_residual_seconds':statistics.mean(r[mode]['rounding_residual_seconds'] for r in selected)}
    return {'steps':rows,'summaries':summaries,
            'interpretation':'An exact accounting of instrumented time, not a causal or recoverable-overhead decomposition. Shared labels do not guarantee matched numerical boundaries. Slowest-rank groups come from the same actual slowest rank per implementation/step; never sum independent per-region maxima.'}


def review(root,ranks=83):
    r=json.loads((root/'results.json').read_text())
    if r['status']!='passed' or not all(g['passed'] for g in r['gates'].values()):
        raise ValueError('Passing correctness campaign required')
    runs=r['runs'];calibration={}
    for kind in ('candidate','legacy'):
        timed=[v['ordinary_step_median'] for v in runs.values() if v['kind']==kind and v['phase']=='timing']
        observed=runs[kind+'-profile']['ordinary_step_median']
        base=statistics.median(timed)
        calibration[kind]={'uninstrumented_median_of_step_medians':base,'profile_step_median':observed,
                           'profile_difference_percent':100*(observed/base-1),'timing_run_count':len(timed)}
    for label in ('candidate-profile','legacy-profile'):
        parsed=analyze(root/label/'run.log')
        if not parsed['completed']:raise ValueError('Incomplete profile run')
        if [s['state'] for s in parsed['step_records']] != [s['state'] for s in analyze(root/'candidate-warmup/run.log')['step_records']]:
            raise ValueError('Different numerical workload')
    b=block_records(root/'candidate-profile/run.log',ranks)
    l=legacy_records(root/'legacy-profile',ranks)
    paired=ledger(b,l)
    wt=compare_profiles(root,ranks)
    whole={k:sum(statistics.mean(x['wall'] for x in rows.values()) for rows in data.values()) for k,data in [('block',b),('legacy',l)]}
    wg=wt['candidate']['wavelet_seconds']-wt['legacy']['wavelet_seconds']
    fullgap=whole['block']-whole['legacy']
    return {'status':'accounted_not_causally_explained','campaign_sha256':hashlib.sha256((root/'results.json').read_bytes()).hexdigest(),
            'headline':{'block_seconds':calibration['candidate']['uninstrumented_median_of_step_medians'],
                        'legacy_seconds':calibration['legacy']['uninstrumented_median_of_step_medians'],
                        'gap_seconds':calibration['candidate']['uninstrumented_median_of_step_medians']-calibration['legacy']['uninstrumented_median_of_step_medians']},
            'profiling_calibration':calibration,'paired_profile_accounting':paired,
            'whole_window_including_checkpoint':{'rank_mean_step_wall_sums':whole,'gap_seconds':fullgap,
                'matched_wavelet_gap_seconds':wg,'outside_wavelet_gap_seconds':fullgap-wg,
                'wavelet_calls':wt['candidate']['wavelet_calls'],
                'scope':'Inclusive matched WT_after_step envelope and all remaining step time. Contains checkpoint; do not divide to infer ordinary-step or headline percentages.'},
            'limitations':['Single detailed profile per implementation; calibration includes run variation and is not an overhead correction.',
                           'Printed CPU timestep medians and rank wall-clock profile means are different statistics; never force them to reconcile.',
                           'Signed profile groups account for totals but do not prove redundant computations or avoidable costs.',
                           'MPI timing can include waiting, progress and local work. Process CPU can include spinning.',
                           'No rank-synchronized event chronology or hardware-counter evidence is present.',
                           'Source/binary provenance and checkpoint reports are audited separately by analyze_stage194_campaign.py.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--ranks',type=int,default=83)
    a=p.parse_args();print(json.dumps(review(a.results,a.ranks),indent=2))
