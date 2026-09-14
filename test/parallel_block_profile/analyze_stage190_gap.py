#!/usr/bin/env python3
"""Attribute ordinary block steps; do not treat block phase time as legacy excess.

Usage: python3 analyze_stage190_gap.py CANDIDATE_PROFILE_LOG RESULTS_JSON
Uses the existing detailed-profile parser; emits a compact, reproducible audit.
"""
import argparse
import collections
import hashlib
import json
import statistics
from pathlib import Path

from analyze import analyze

GROUPS = {
    'Domain boundary handling': [28, 38, 58, 59, 60, 61, 62, 63],
    'Block reconstruction and migration': [5, 6, 7, 8, 9, 10, 11, 12, 13, 35, 48],
    'Native inverse': [14, 15, 16, 17, 18, 19, 20, 21, 36, 41, 42],
    'Scalar input and restriction': [22, 23, 24, 25, 37, 51],
    'Block ghost handling': [33, 34],
    'Native mass and velocity': [30, 31, 39],
    'Shared Domain producers': [26, 27, 29, 52],
    'Domain/block publication': [32, 40, 49, 50],
    'Wavelet and compression': [43, 44, 65],
    'RK and tendency assembly': [45, 46, 47, 53],
    'Other timestep work': [1, 2, 3, 4, 54, 55, 56, 57, 64],
}
MPI_IDS = [13, 19, 34, 37, 38, 39, 40]


def audit(log, results):
    parsed = analyze(log)
    campaign = json.loads(results.read_text())
    if not parsed['completed'] or campaign['status'] != 'passed':
        raise ValueError('Completed profile and passing campaign required')
    ids = [i for group in GROUPS.values() for i in group]
    assert sorted(ids) == list(range(1, 66))
    windows = parsed['detail_windows']
    if any(w.get('region_count') != 65 for w in windows):
        raise ValueError('Expected the Stage 190 65-region schema')
    steps = [s for w in windows for s in w['step_summaries']]
    if not steps or any(not s['complete'] for s in steps):
        raise ValueError('Missing or dropped rank/step records')
    ordinary = [s for s in steps if not s['restart']]
    if not ordinary:
        raise ValueError('No ordinary steps')
    if any(s['self_avg'][56] != 0 for s in steps):
        raise ValueError('Oracle time is present; this is not the optimized production profile')
    if max(s['max_self_conservation_error'] for s in ordinary) > 1e-4:
        raise ValueError('Exclusive timings do not conserve timestep time')
    groups = {}
    for name, group in GROUPS.items():
        groups[name] = {
            'region_ids': group,
            'mean_rank_mean_seconds': statistics.mean(
                sum(s['self_avg'][i-1] for i in group) for s in ordinary),
            'mean_seconds_on_each_steps_slowest_rank': statistics.mean(
                sum(s['self_on_max_wall_rank'][i-1] for i in group) for s in ordinary),
        }
    rank_records = collections.defaultdict(list)
    for w in windows:
        for s in w['steps'].values():
            if s['restarts']:
                continue
            mpi = sum(s['self-wall'][i-1] for i in MPI_IDS)
            rank_records[s['rank']].append((s['wall'], mpi, s['wall']-mpi))
    if len(rank_records) != campaign['ranks'] or any(len(v) != len(ordinary) for v in rank_records.values()):
        raise ValueError('Rank coverage differs from campaign')
    ranks = {rank: dict(zip(('mean_step_seconds', 'mean_named_MPI_seconds',
                            'mean_outside_named_MPI_seconds'),
                           (statistics.mean(v[i] for v in records) for i in range(3))))
             for rank, records in sorted(rank_records.items())}
    counters = collections.Counter()
    boundary = collections.defaultdict(lambda: {'wall_rank_mean_sum': 0.0, 'global_calls': 0})
    names = {r['id']: r['region'] for w in windows for r in w['regions']}
    for w in windows:
        counters.update({k.strip(): v for k, v in w['counters'].items()})
        for caller in w['boundary_callers']:
            dest = boundary[names.get(caller['caller_id'], str(caller['caller_id']))]
            dest['wall_rank_mean_sum'] += caller['wall_avg']
            dest['global_calls'] += caller['calls']
    return {
        'status': 'passed_profile_audit',
        'source_commit': 'cdd2eb64',
        'inputs_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (log, results)},
        'campaign_summary': campaign['summary'],
        'ordinary_steps': [s['sequence'] for s in ordinary],
        'ordinary_steps_with_remap': [s['sequence'] for s in ordinary if s['remap']],
        'profile_mean_rank_mean_step_seconds': statistics.mean(s['wall_avg'] for s in ordinary),
        'profile_mean_slowest_rank_step_seconds': statistics.mean(s['wall_max'] for s in ordinary),
        'max_step_conservation_error': max(s['max_self_conservation_error'] for s in ordinary),
        'groups': groups,
        'proof_only_mean_seconds': statistics.mean(sum(s['self_avg'][i-1] for i in (6, 7, 9)) for s in ordinary),
        'mean_named_MPI_seconds': statistics.mean(v['mean_named_MPI_seconds'] for v in ranks.values()),
        'named_MPI_region_ids': MPI_IDS,
        'ranks': ranks,
        'whole_profile_window_counters': dict(counters),
        'whole_profile_window_boundary_callers': dict(boundary),
        'limitations': [
            'The campaign ratio uses unprofiled medians; phase attribution uses means of 11 profiled non-checkpoint steps.',
            'Remap steps remain ordinary, matching the campaign definition.',
            'Named MPI times include spin and arrival imbalance; outside-named-MPI time may include unclassified MPI.',
            'Phase groups are disjoint block costs, not measured excess over the legacy phases.',
            'No legacy phase profile or controlled ablation was supplied. Do not claim a causal decomposition of the entire gap.',
            'Whole-window counters and boundary caller totals include the checkpoint step and possible between-step work.',
        ],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate_profile', type=Path)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.candidate_profile, args.results), indent=2))
