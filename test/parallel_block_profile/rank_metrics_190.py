#!/usr/bin/env python3
"""Linux Slurm placement and per-rank peak RSS; no polling during timed solver work."""
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys


def placement():
    return {'memory_allowed':next((l.split(':',1)[1].strip() for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('Mems_allowed_list:')),None),
            'rank':int(os.environ.get('SLURM_PROCID',0)), 'host':socket.gethostname(),
            'allowed_cpus':sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None}


def memory():
    result=placement()
    for name in ('meminfo','vmstat'):
        text=Path('/proc',name).read_text()
        result[name]={parts[0].rstrip(':'):int(parts[1]) for line in text.splitlines()
                      if len(parts:=line.split())>=2 and parts[1].isdigit()}
    return result


def collect(mode,ranks):
    command=['srun',f'--ntasks={ranks}','--cpus-per-task=1','--cpu-bind=cores','--kill-on-bad-exit=1']
    if mode=='--memory': command+=['--ntasks-per-node=1']
    command += [sys.executable,str(Path(__file__).resolve()),mode]
    run=subprocess.run(command,capture_output=True,text=True,check=True)
    rows=[json.loads(line) for line in run.stdout.splitlines() if line.startswith('{')]
    if len(rows)!=ranks: raise ValueError('Missing rank/node placement measurements')
    return rows


def snapshot_allocation(out,ranks):
    result={'ranks':collect('--probe',ranks),'scope':'Actual affinity under the same srun core binding as measured runs.'}
    placements=result['ranks']
    occupied=set()
    for row in placements:
        cpus=row['allowed_cpus']
        if not cpus: raise ValueError('No CPU affinity reported')
        for cpu in cpus:
            key=(row['host'],cpu)
            if key in occupied: raise ValueError('Overlapping rank CPU affinities')
            occupied.add(key)
    topology=subprocess.run(['lscpu','--parse=CPU,CORE,SOCKET,NODE'],capture_output=True,text=True)
    result['controller_cpu_topology']=topology.stdout
    run=subprocess.run(['lscpu'],capture_output=True,text=True)
    result['controller_lscpu']=run.stdout
    result['controller_hostname']=socket.gethostname()
    (out/'placement.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def node_memory(path):
    rows=collect('--memory',int(os.environ.get('SLURM_JOB_NUM_NODES','1')))
    path.write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def rank_resources(directory):
    rows=[]
    for path in sorted(directory.glob('rank-*.json')):
        row=json.loads(path.read_text())
        text=path.with_suffix('.resources').read_text()
        for key,phrase in [('max_rss_kib','Maximum resident set size (kbytes)'),
                           ('major_faults','Major (requiring I/O) page faults'),
                           ('minor_faults','Minor (reclaiming a frame) page faults')]:
            match=re.search(re.escape(phrase)+r':\s*(\d+)',text)
            if not match: raise ValueError('Missing child resource count: '+str(path))
            row[key]=int(match[1])
        rows.append(row)
    return rows


if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='--probe': print(json.dumps(placement()),flush=True)
    elif mode=='--memory': print(json.dumps(memory()),flush=True)
    else: raise SystemExit('Unknown rank metrics mode')
