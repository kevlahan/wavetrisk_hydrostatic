"""Explicit local J5 numerical screening policy, not byte-equality or a proof."""
from collections import Counter, deque
import re
import numpy as np
from compare_restart_checkpoints import read, compare

# Fixed before the lighter-fixture experiments; do not auto-calibrate to a
# candidate. Existing scientific oracle tolerances are untouched. These
# absolute limits are a conservative screen around the independently validated
# 179d double-precision checkpoint differences, not a universal accuracy bound.
POLICY={'id':'stage180-local-v1','atmosphere_absolute':{
    'velocity':1e-10,'mass':1e-10,'temperature':1e-7},'soil_absolute':0.0}


def topology(path):
    counts=Counter();nonzero=0
    for _,threshold,poles,coarse,waves,required in read(path,-10,30):
        if not all(np.isfinite(a).all() for a in (threshold,poles,coarse,waves)):
            raise ValueError('Nonfinite fixture/checkpoint')
        queue=deque([5]*4)
        for flags in required:
            if not queue: raise ValueError('Trailing checkpoint wavelet records')
            level=queue.popleft();counts[level]+=1
            if any(x not in (0,1) for x in flags): raise ValueError('Unsupported logical representation')
            queue.extend([level+1]*int(np.count_nonzero(flags)))
        if queue: raise ValueError('Missing checkpoint wavelet records')
        nonzero+=int(np.count_nonzero(waves))
    if not counts or not nonzero: raise ValueError('Empty or unevolved checkpoint')
    return {'patches_by_level':dict(sorted(counts.items())),'maximum_level':max(counts),
            'nonzero_wavelet_values':nonzero}


def accept(report,exact=False):
    reasons=[]
    for key in ('header_differences','topology_differences','threshold_nonfinite'):
        if report.get(key,0): reasons.append(key)
    if report['threshold_max_abs']!=0: reasons.append('threshold values differ')
    for name,field in report['fields'].items():
        _,variable,zone=name.split('/')
        limit=0.0 if exact or zone=='soil' else POLICY['atmosphere_absolute'][variable]
        if field['nonfinite']: reasons.append(name+': nonfinite')
        if not np.isfinite(field['max_abs']) or field['max_abs']>limit:
            reasons.append(f"{name}: {field['max_abs']:.17g} exceeds {limit:.17g}")
    if not report['fields']: reasons.append('No compared fields')
    return reasons


def compare_runs(reference,candidate,seed_name,exact=False):
    a={p.name:p for p in reference.glob('*_checkpoint_*.bin.zst') if p.name!=seed_name}
    b={p.name:p for p in candidate.glob('*_checkpoint_*.bin.zst') if p.name!=seed_name}
    if a.keys()!=b.keys() or len(a)<2:
        raise ValueError('Require matching sets of at least two newly written checkpoints')
    failures=[];reports={}
    for name in sorted(a):
        left=topology(a[name]);right=topology(b[name])
        if left['patches_by_level']!=right['patches_by_level']: failures.append(name+': patch levels differ')
        report=compare(a[name],b[name]);reports[name]=report
        failures.extend(name+': '+r for r in accept(report,exact))
    for root in (reference,candidate):
        log=(root/'run.log').read_text()
        if 'Remapping vertical coordinates' not in log: failures.append(str(root)+': no remap event')
        for name in a:
            index=int(name.split('_checkpoint_')[1].split('.')[0])
            if not re.search(r'Restarting from checkpoint\s+'+str(index)+r'\b',log):
                failures.append(str(root)+': missing reload '+str(index))
    return {'passed':not failures,'failures':failures,'policy':POLICY,'exact_fields_required':exact,'checkpoints':reports}
