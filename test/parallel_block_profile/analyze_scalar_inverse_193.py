#!/usr/bin/env python3
"""Validate scalar-inverse work counts and compare matched transform profiles."""
import re
from pathlib import Path
from analyze import analyze
from analyze_inverse_192 import compare_profiles as compare_wavelet


def scalar_work(path):
    rows=re.findall(r'^\s*scalar inverse columns: stencils samples addresses geometry old-addresses old-geometry = ([0-9 ]+)$',
                    Path(path).read_text(),re.M)
    values=[list(map(int,r.split())) for r in rows]
    if not values or any(len(v)!=6 for v in values): raise ValueError('Missing scalar inverse work counters')
    totals=[sum(v[k] for v in values) for k in range(6)]
    names=('compiled_stencils','produced_field_samples','resolved_field_addresses','geometry_reads',
           'reference_equivalent_field_addresses','reference_equivalent_geometry_reads')
    report=dict(zip(names,totals))
    if not all(totals) or totals[4]<=totals[2] or totals[5]<=totals[3]:
        raise ValueError('Scalar inverse columns did not reduce preparation work')
    report['address_reduction_percent']=100*(1-totals[2]/totals[4])
    report['geometry_read_reduction_percent']=100*(1-totals[3]/totals[5])
    report['scope']='Global sums over all profile windows, including checkpoint. Reference-equivalent counts derive from the retained scalar kernel; arithmetic sample count is unchanged. Field-address counts exclude the separate geometry accessor lookups.'
    return report


def compare_profiles(root,ranks=83):
    result=compare_wavelet(root,ranks)
    result['scalar_work']=scalar_work(Path(root)/'candidate-profile/run.log')
    # A computational-only experiment must retain the same transfer schedule.
    for family,b in result['baseline']['transfers'].items():
        c=result['candidate']['transfers'][family]
        for key in ('remote_bytes','remote_messages','calls_per_rank'):
            if b[key]!=c[key]: raise ValueError('Transfer workload changed: '+family+' '+key)
    kernels={}
    for kind in ('baseline','candidate'):
        rows=[r for w in analyze(Path(root)/(kind+'-profile/run.log'))['detail_windows']
              for r in w['regions'] if r['id']==36]
        if not rows: raise ValueError('Missing inverse kernel profile')
        kernels[kind]={'inclusive_seconds':sum(r['inclusive_wall_avg'] for r in rows),
                       'global_calls':sum(r['global_calls'] for r in rows)}
    if kernels['baseline']['global_calls']!=kernels['candidate']['global_calls']:
        raise ValueError('Inverse kernel call coverage changed')
    kernels['scope']='Matched scalar plus inner-vector kernel scopes, whole-window rank-mean wall time. Includes address/geometry preparation and arithmetic, excludes inverse MPI. Not added to parent WT or inverse time.'
    result['inverse_local_kernels']=kernels
    result['transfer_workload_unchanged']=True
    return result
