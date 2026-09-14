#!/usr/bin/env python3
"""Compare complete WT_after_step scopes, not disjoint inverse/boundary labels."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import statistics
from analyze import analyze
from analyze_stage191_campaign import audit, TESTED_SOURCE


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def block_regions(path):
    a=analyze(path);regions={}
    for window in a['detail_windows']:
        for item in window['regions']:
            r=regions.setdefault(item['id'],{'name':item['region'],'self_wall':0.,'inclusive_wall':0.,'self_cpu':0.,'calls':0})
            for key,field in [('self_wall','self_wall_avg'),('inclusive_wall','inclusive_wall_avg'),('self_cpu','self_cpu_avg'),('calls','global_calls')]:r[key]+=item[field]
    return regions,a


def legacy_regions(root,ranks):
    regions={};calls=defaultdict(list)
    paths=sorted(root.glob('legacy-detail-rank-*.txt'));assert len(paths)==ranks
    for path in paths:
        for line in path.read_text().splitlines():
            f=line.split()
            if not f or f[0]!='region':continue
            i=int(f[1]);r=regions.setdefault(i,{'self_wall':0.,'inclusive_wall':0.,'self_cpu':0.,'calls':0})
            r['self_wall']+=float(f[2])/ranks;r['self_cpu']+=float(f[3])/ranks
            r['inclusive_wall']+=float(f[4])/ranks;r['calls']+=int(f[6]);calls[i].append(int(f[6]))
    assert set(calls[65])=={60} and set(calls[1])=={12}
    return regions


def native_transfers(path):
    records={}
    for line in path.read_text().splitlines():
        if not line.strip().startswith(('inverse native gather:','inverse native scatter:','inverse native aliases:')):continue
        f=line.split('|');name,count=f[0].strip().split(':')
        assert len(f)==6
        r=records.setdefault(name,{'rank_mean_seconds':0.,'reported_rank_mean_calls':0,'remote_messages':0,'remote_payload_bytes':0,'windows':0})
        r['rank_mean_seconds']+=float(f[1].split()[1]);r['reported_rank_mean_calls']+=int(count)
        r['remote_messages']+=int(f[4]);r['remote_payload_bytes']+=int(f[5]);r['windows']+=1
    assert len(records)==3 and all(v['windows']==2 for v in records.values())
    return records


def compare(root,repo):
    checked=audit(root,repo)
    results=root/'results';ranks=83;steps=12
    candidate,parsed=block_regions(results/'candidate-profile/run.log')
    baseline,_=block_regions(results/'baseline-profile/run.log')
    legacy=legacy_regions(results/'legacy-profile',ranks)
    assert candidate[65]['calls']==baseline[65]['calls']==legacy[65]['calls']==60*ranks
    assert candidate[1]['calls']==baseline[1]['calls']==legacy[1]['calls']==steps*ranks
    wt={name:{'whole_window_rank_mean_seconds':v[65]['inclusive_wall'],
              'seconds_per_step_including_checkpoint':v[65]['inclusive_wall']/steps,
              'calls_per_rank':v[65]['calls']/ranks}
        for name,v in [('stage190',baseline),('stage191',candidate),('legacy',legacy)]}
    excess=candidate[65]['inclusive_wall']-legacy[65]['inclusive_wall']
    step_excess=candidate[1]['inclusive_wall']-legacy[1]['inclusive_wall']
    transfers=native_transfers(results/'candidate-profile/run.log')
    total={k:sum(v[k] for v in transfers.values()) for k in ('rank_mean_seconds','remote_messages','remote_payload_bytes','reported_rank_mean_calls')}
    assert total['reported_rank_mean_calls']*ranks==candidate[19]['calls']
    assert total['rank_mean_seconds']<candidate[14]['inclusive_wall']<candidate[65]['inclusive_wall']
    rank_data=defaultdict(list);mpi_ids=[13,19,34,37,38,39,40]
    for w in parsed['detail_windows']:
        for s in w['steps'].values():
            if s['restarts']:continue
            mpi=sum(s['self-wall'][i-1] for i in mpi_ids)
            rank_data[s['rank']].append((s['wall'],mpi,s['wall']-mpi))
    assert len(rank_data)==ranks and all(len(v)==11 for v in rank_data.values())
    ranks_summary={k:dict(zip(('wall','named_MPI','outside_named_MPI'),[statistics.mean(row[i] for row in rows) for i in range(3)])) for k,rows in sorted(rank_data.items())}
    return {'status':'passed_matched_scope_audit','job':54561,'tested_source_commit':TESTED_SOURCE,
        'input_campaign_sha256':checked['results_sha256'],
        'scope':'Whole-window inclusive WT_after_step wall time, averaged across ranks, 60 calls/rank and 12 timesteps including checkpoint. No sum of phase maxima.',
        'matched_wavelet_transaction':wt,'wavelet_ratio':candidate[65]['inclusive_wall']/legacy[65]['inclusive_wall'],
        'wavelet_excess_whole_window_seconds':excess,'wavelet_excess_seconds_per_step_including_checkpoint':excess/steps,
        'whole_profile_step_excess_seconds':step_excess,'fraction_of_whole_profile_step_excess_in_wavelet_percent':100*excess/step_excess,
        'native_inverse':{'calls_per_rank':candidate[14]['calls']/ranks,'inclusive_rank_mean_seconds':candidate[14]['inclusive_wall'],
            'explicit_wait_self_seconds':candidate[19]['self_wall'],'explicit_wait_process_cpu_seconds':candidate[19]['self_cpu'],
            'regular_plus_outer_kernel_inclusive_seconds':candidate[36]['inclusive_wall'],
            'transport_helpers':transfers,'transport_helpers_total':total,
            'transport_helper_fraction_of_inverse_percent':100*total['rank_mean_seconds']/candidate[14]['inclusive_wall']},
        'ordinary_rank_diagnostic':ranks_summary,
        'source_findings':['WT_after_step contains the entire forward/compression/inverse transaction in both implementations; inclusive timing already charges nested legacy boundary work.',
        'Native output skips legacy forward-wavelet loops and legacy inverse; repeated whole-transform call counts are equal. No evidence of a duplicate full legacy transform on each native stage.',
        'Regular inverse kernels stay on final block owners; outer-edge dependencies retain geometric route owners. transfer_nodes gathers/scatters between these layouts, and aliases use a separate exchange.',
        'The provisional inverse also gathers/scatters boundary state and publishes Domain consumer state after the interior computation.',
        'Native gather/scatter/alias payload counters count remote sends only; they exclude self copies and do not represent total DRAM traffic.'],
        'next_experiment':'Within the inverse transaction, remove or narrow per-level final-owner/geometric-owner transfers while preserving outer-edge write order and consumer aliases; measure the same full transaction against legacy.',
        'remaining_limits':['This isolates a major cost center, not the whole 0.783 s unprofiled ordinary-step excess.',
        'Whole-window inclusive records cannot exclude the checkpoint step or give per-step parent/child attribution. Do not divide this excess by the unprofiled ordinary-step excess.',
        'Transport helpers include packing, copying and waiting; their sum is contained in inverse inclusive time and must not be added to it.',
        'MPI process CPU includes spin; outside-named-MPI time is not guaranteed useful arithmetic. Arrival imbalance and implementation costs both remain relevant.',
        'Legacy inverse timer has 73 calls/rank versus 48 native inverse calls, so their bare timer rows are not matched transactions. The common 60-call WT parent is the defensible comparison.',
        'Proof of recoverable savings requires a controlled ownership-transfer ablation; payload volume alone is not that proof.']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('extracted_archive',type=Path);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
    print(json.dumps(compare(a.extracted_archive,a.repo),indent=2))
