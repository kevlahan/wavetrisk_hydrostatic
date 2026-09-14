#!/usr/bin/env python3
"""Compare matched complete wavelet transactions and inverse transfer traffic."""
import argparse
import json
import math
from pathlib import Path
from analyze import analyze
from legacy_metrics_191 import collect


def block_profile(path, ranks):
    parsed=analyze(path)
    if not parsed['completed']: raise ValueError('Incomplete block profile')
    regions={}
    for window in parsed['detail_windows']:
        for item in window['regions']:
            row=regions.setdefault(item['id'],{'inclusive_wall':0.,'calls':0})
            row['inclusive_wall']+=item['inclusive_wall_avg']
            row['calls']+=item['global_calls']
    if not all(i in regions for i in (1,14,19,65)): raise ValueError('Missing matched profile regions')
    steps=len(parsed['step_records'])
    if regions[1]['calls']!=steps*ranks: raise ValueError('Incomplete block step/rank coverage')
    traffic={}
    for line in path.read_text().splitlines():
        if not line.strip().startswith(('inverse native gather:','inverse native scatter:','inverse native aliases:')): continue
        fields=line.split('|')
        if len(fields)!=6: raise ValueError('Unsupported transfer record')
        name,calls=fields[0].strip().split(':')
        row=traffic.setdefault(name,{'seconds':0.,'calls_per_rank':0,'remote_messages':0,'remote_bytes':0,'windows':0})
        row['seconds']+=float(fields[1].split()[1]);row['calls_per_rank']+=int(calls)
        row['remote_messages']+=int(fields[4]);row['remote_bytes']+=int(fields[5]);row['windows']+=1
    if len(traffic)!=3 or any(v['windows']!=len(parsed['detail_windows']) for v in traffic.values()):
        raise ValueError('Missing transfer profile window')
    if sum(v['calls_per_rank'] for v in traffic.values())*ranks!=regions[19]['calls']:
        raise ValueError('Transfer/MPI completion coverage differs')
    for row in traffic.values():
        if any(not math.isfinite(v) or v<0 for v in row.values()): raise ValueError('Invalid transfer metric')
        row['remote_bytes_per_step_including_checkpoint']=row['remote_bytes']/steps
        row['remote_messages_per_step_including_checkpoint']=row['remote_messages']/steps
    if any(not math.isfinite(v['inclusive_wall']) or v['inclusive_wall']<0 for v in regions.values()):
        raise ValueError('Invalid inclusive time')
    return {'steps':steps,'wavelet_calls':regions[65]['calls'],
            'wavelet_seconds':regions[65]['inclusive_wall'],'inverse_seconds':regions[14]['inclusive_wall'],
            'transfers':traffic}


def compare_profiles(root,ranks=83):
    root=Path(root)
    baseline=block_profile(root/'baseline-profile/run.log',ranks)
    candidate=block_profile(root/'candidate-profile/run.log',ranks)
    legacy=collect(root/'legacy-profile',ranks)
    calls=[];seconds=0.
    for path in sorted((root/'legacy-profile').glob('legacy-detail-rank-*.txt')):
        rows=[line.split() for line in path.read_text().splitlines() if line.startswith('region')]
        wt=[f for f in rows if int(f[1])==65]
        if len(wt)!=1: raise ValueError('Missing/duplicate legacy WT scope')
        calls.append(int(wt[0][6]));seconds+=float(wt[0][4])/ranks
    if (baseline['steps']!=candidate['steps'] or baseline['steps']!=len(legacy['steps']) or
        baseline['wavelet_calls']!=candidate['wavelet_calls'] or len(set(calls))!=1 or
        sum(calls)!=candidate['wavelet_calls'] or not math.isfinite(seconds) or seconds<=0):
        raise ValueError('Unmatched complete wavelet transaction')
    b=baseline['transfers']['inverse native scatter']['remote_bytes']
    c=candidate['transfers']['inverse native scatter']['remote_bytes']
    return {'baseline':baseline,'candidate':candidate,
            'legacy':{'wavelet_seconds':seconds,'wavelet_calls':sum(calls),'steps':len(legacy['steps'])},
            'scatter_byte_reduction_percent':100*(1-c/b) if b else None,
            'wavelet_reduction_percent':100*(1-candidate['wavelet_seconds']/baseline['wavelet_seconds']),
            'candidate_wavelet_to_legacy_ratio':candidate['wavelet_seconds']/seconds,
            'scope':'Whole-window rank-mean inclusive WT_after_step, including checkpoint. Helpers are contained in inverse time. Remote sends exclude self copies; bytes are not measured DRAM or NIC traffic. Profiles are excluded from headline timing. Per-step traffic here includes checkpoint; it is not an ordinary-step measurement.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--ranks',type=int,default=83)
    a=p.parse_args();print(json.dumps(compare_profiles(a.results,a.ranks),indent=2))
