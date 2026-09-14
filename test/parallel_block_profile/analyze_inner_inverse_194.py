#!/usr/bin/env python3
"""Separate scalar, inner and outer inverse work without adding nested timers."""
import math
from pathlib import Path
import re
import statistics
from analyze import analyze
from analyze_scalar_inverse_193 import compare_profiles as compare_previous


def inner_work(path):
    text=Path(path).read_text()
    prefix='inner inverse geometry: prepared-cells cell-updates reads reference-reads = '
    rows=re.findall(r'^\s*'+re.escape(prefix)+r'([0-9 ]+)$',text,re.M)
    if not rows or len(rows)!=text.count(prefix): raise ValueError('Missing/malformed inner geometry counters')
    rows=[list(map(int,r.split())) for r in rows]
    if any(len(r)!=4 for r in rows): raise ValueError('Invalid inner geometry counter width')
    totals=[sum(r[i] for r in rows) for i in range(4)]
    prepared,cells,reads,reference=totals
    if cells<=0 or reference!=32*cells or reads<=0 or reads>reference:
        raise ValueError('Invalid inner geometry work coverage')
    if prepared and reads!=32*prepared: raise ValueError('Prepared cells/geometry reads disagree')
    if not prepared and reads!=reference: raise ValueError('Uncached geometry coverage differs')
    return dict(zip(('prepared_cells','physical_cell_updates','geometry_reads','reference_geometry_reads'),totals),
                geometry_read_reduction_percent=100*(1-reads/reference))


def split_profile(path):
    a=analyze(Path(path))
    if not a['completed'] or not a['detail_windows']: raise ValueError('Incomplete split profile')
    phases={21:'outer_vector',36:'combined_parent',66:'scalar',67:'inner_vector'}
    total={name:{'inclusive_wall_seconds':0.,'self_wall_seconds':0.,'self_cpu_seconds':0.,'global_calls':0}
           for name in phases.values()}
    for w in a['detail_windows']:
        rows={r['id']:r for r in w['regions']}
        if not set(phases)<=set(rows): raise ValueError('Missing scalar/inner/outer split')
        for i,name in phases.items():
            r=rows[i];v=total[name]
            for key,source in [('inclusive_wall_seconds','inclusive_wall_avg'),('self_wall_seconds','self_wall_avg'),
                               ('self_cpu_seconds','self_cpu_avg'),('global_calls','global_calls')]:
                if not math.isfinite(r[source]) or r[source]<0: raise ValueError('Invalid split timing')
                v[key]+=r[source]
        children=sum(rows[i]['inclusive_wall_avg'] for i in (21,66,67))
        if abs(rows[36]['inclusive_wall_avg']-rows[36]['self_wall_avg']-children)>1e-4:
            raise ValueError('Inverse child timings do not conserve the parent')
        if rows[36]['global_calls']!=sum(rows[i]['global_calls'] for i in (21,66,67)):
            raise ValueError('Inverse split call coverage differs')
        if rows[66]['global_calls']!=2*rows[67]['global_calls']:
            raise ValueError('Scalar phase/inner call coverage differs')
    steps=[s for w in a['detail_windows'] for s in w['step_summaries'] if not s.get('restart')]
    if not steps or any(not s['complete'] or s['max_self_conservation_error']>1e-4 for s in steps):
        raise ValueError('Incomplete ordinary-step split')
    ordinary={name:statistics.mean(s['self_avg'][i-1] for s in steps) for i,name in phases.items()}
    ordinary.pop('combined_parent')
    ordinary['parent_residual']=statistics.mean(s['self_avg'][35] for s in steps)
    return {'whole_window_rank_mean':total,'ordinary_step_mean_rank_mean_self_seconds':ordinary,
            'ordinary_steps':len(steps),
            'scope':'Whole-window times include checkpoint. Parent includes scalar, inner and outer; do not add inclusive parent to its children. Ordinary-step means are instrumented diagnostics, not headline timings.'}


def compare_profiles(root,ranks=83):
    root=Path(root)
    result=compare_previous(root,ranks)
    # The Stage 193 executable retains its original profile schema. The same
    # Stage 194 executable with geometry disabled provides the split control.
    control=split_profile(root/'candidate-control-profile/run.log')
    candidate=split_profile(root/'candidate-profile/run.log')
    work_control=inner_work(root/'candidate-control-profile/run.log')
    work=inner_work(root/'candidate-profile/run.log')
    if work_control['prepared_cells'] or not work['prepared_cells']:
        raise ValueError('Wrong geometry control/candidate selector')
    if work['physical_cell_updates']!=work_control['physical_cell_updates'] or work['reference_geometry_reads']!=work_control['geometry_reads']:
        raise ValueError('Different inner-vector work in control and candidate')
    for name in ('combined_parent','scalar','inner_vector','outer_vector'):
        if control['whole_window_rank_mean'][name]['global_calls']!=candidate['whole_window_rank_mean'][name]['global_calls']:
            raise ValueError('Split comparison has unmatched call counts')
    result['inner_geometry_comparison']={'control':control,'candidate':candidate,'control_work':work_control,'candidate_work':work}
    result['inverse_local_kernels']['scope']='All inverse local kernels, including the outer-vector child; historical combined timer ID 36.'
    return result
