#!/usr/bin/env python3
"""Summarize gated local results without summing unrelated rank maxima."""
import argparse
import json
from pathlib import Path
import statistics
from local_experiment import legacy_detail
from local_campaign import ordinary


def aggregate(steps,names):
    if not steps: raise ValueError('No ordinary complete profile steps')
    means=[statistics.mean(s['self'][i] for s in steps) for i in range(65)]
    return {'ordinary_steps':len(steps),
            'mean_max_rank_wall':statistics.mean(s['wall'] for s in steps),
            'mean_cpu_on_max_wall_rank':statistics.mean(s['cpu'] for s in steps),
            'attribution':'Exclusive wall on each timestep\'s actual slowest rank; not independent region maxima',
            'regions':sorted([{'id':i+1,'name':names.get(i+1,str(i+1)),
                               'self_wall_on_step_critical_rank':value} for i,value in enumerate(means)],
                             key=lambda r:r['self_wall_on_step_critical_rank'],reverse=True)}


def profiles(result):
    record=result['validation']['block-detail'];steps=[];names={}
    for window in record['analysis']['detail_windows']:
        names.update({r['id']:r['region'] for r in window['regions']})
        for summary in window['step_summaries']:
            if not summary['complete']: raise ValueError('Incomplete block profile window')
            if summary['max_self_conservation_error']>1e-6: raise ValueError('Profile conservation failed')
            if summary['restart'] or summary['remap']: continue
            rank=summary['max_wall_rank']
            local=window['steps'][str(summary['sequence'])+':'+str(rank)]
            steps.append({'wall':summary['wall_max'],'cpu':local['cpu'],
                          'self':summary['self_on_max_wall_rank']})
    native=aggregate(steps,names)
    directory=Path(result['validation']['legacy-detail']['directory'])
    ranks=legacy_detail(directory,result['binding']['ranks'])
    sequences=set(ranks[0]['steps'])
    if any(set(r['steps'])!=sequences for r in ranks): raise ValueError('Legacy rank windows differ')
    steps=[]
    for sequence in sorted(sequences):
        rows=[r['steps'][sequence] for r in ranks]
        if any(r['restart'] or r['remap'] for r in rows):continue
        slow=max(rows,key=lambda r:r['wall'])
        steps.append({'wall':slow['wall'],'cpu':slow['cpu'],'self':[slow['self'][i][0] for i in range(1,66)]})
    return {'block':native,'legacy':aggregate(steps,names),
            'caution':'Shared region IDs do not guarantee identical implementation boundaries; no phase subtraction or speedup extrapolation.'}


def gate_summary(gate):
    maxima={}
    for report in gate['checkpoints'].values():
        for name,field in report['fields'].items():
            maxima[name]=max(maxima.get(name,0),field['max_abs'])
    return {'passed':gate['passed'],'failures':gate['failures'],'max_abs_by_field':maxima}


def summarize(path):
    result=json.loads(path.read_text())
    summary={'scope':result['scope'],'binding':result['binding'],'validation':{},'timing':[]}
    for label,record in result['validation'].items():
        summary['validation'][label]={**gate_summary(record['gate']),
                                     'memory':{k:v for k,v in record['memory'].items() if k!='samples'}}
    for record in result['timing']:
        summary['timing'].append({'pair':record['pair'],'label':record['label'],
                                  'ordinary_steps':ordinary(record),
                                  'ordinary_mean':statistics.mean(ordinary(record)),
                                  'memory':{k:v for k,v in record['memory'].items() if k!='samples'}})
    summary['paired_ratios']=result.get('paired_ratios',[])
    summary['timing_blocked_reasons']=result.get('timing_blocked_reasons',{})
    for name in ('oracle-diagnostic','legacy-instrumentation-diagnostic'):
        diagnostic=path.parent/(name+'.json')
        if diagnostic.exists():
            summary[name]={label:gate_summary(gate) for label,gate in json.loads(diagnostic.read_text()).items()}
    if (all(g['passed'] for g in summary['validation'].values()) and
            {'block-detail','legacy-detail'}.issubset(result['validation'])):
        summary['profiles']=profiles(result)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    print(json.dumps(summarize(p.parse_args().results),indent=2))
