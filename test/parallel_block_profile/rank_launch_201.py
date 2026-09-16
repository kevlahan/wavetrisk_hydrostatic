#!/usr/bin/env python3
"""Use Stage200 affinity; optionally sample a frozen rank cohort on unchanged optimized observers."""
import json,os,subprocess,sys
from pathlib import Path
from rank_launch_190 import execute,placement

def perf_command(path,command):
    return ['perf','record','--clockid','mono','--no-inherit','-e','cycles:u','-e','instructions:u',
            '-c','20000000','--call-graph','dwarf,8192','-o',str(path),'--',*command]

def export_perf(path):
    with path.with_suffix('.txt').open('x') as out,path.with_suffix('.export.log').open('x') as err:
        subprocess.run(['perf','script','--ns','--show-lost-events','-i',str(path),'-F','pid,time,event,ip,sym,dso'],stdout=out,stderr=err,check=True)
    with path.with_suffix('.header.txt').open('x') as out:
        subprocess.run(['perf','report','--header-only','--stdio','-i',str(path)],stdout=out,stderr=subprocess.STDOUT,check=True)

def main():
    case=json.loads(Path(sys.argv[1]).read_text());rank=int(os.environ['SLURM_PROCID']);cpu=case['mapping'][rank]['cpu']
    if cpu not in os.sched_getaffinity(0):raise ValueError('Requested CPU unavailable')
    os.sched_setaffinity(0,{cpu});os.environ['WAVETRISK_EXPECT_CPU200']=str(cpu);os.environ['MV2_ENABLE_AFFINITY']='0'
    if sys.argv[2]=='--probe':print(json.dumps(placement(rank,[])),flush=True);return 0
    directory=Path(sys.argv[2]);command=['numactl','--localalloc',*sys.argv[3:]]
    sample=rank in json.loads(os.environ.get('WAVETRISK_SAMPLE_RANKS201','[]'))
    if not sample:return execute(directory,command)
    path=directory/f'perf-rank-{rank}.data'
    code=execute(directory,perf_command(path,command))
    resource=directory/f'rank-{rank}.json';r=json.loads(resource.read_text())
    r['resource_scope']='perf wrapper and descendants; not comparable to bare solver-only RSS/CPU metrics'
    resource.write_text(json.dumps(r,indent=2)+'\n')
    if code:return code
    export_perf(path)
    return 0
if __name__=='__main__':raise SystemExit(main())
