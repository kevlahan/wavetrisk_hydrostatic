#!/usr/bin/env python3
"""Pin before MPI initialization; local NUMA allocation; preserve solver wait4 metrics."""
import json,os,sys
from pathlib import Path
from rank_launch_190 import execute,placement

def main():
    case=json.loads(Path(sys.argv[1]).read_text());rank=int(os.environ['SLURM_PROCID']);row=case['mapping'][rank]
    cpu=row['cpu'];allowed=os.sched_getaffinity(0)
    if cpu not in allowed:raise ValueError(f'OS CPU {cpu} unavailable to rank {rank}; step cgroup/affinity restricts the allocation')
    os.sched_setaffinity(0,{cpu})
    if os.sched_getaffinity(0)!={cpu}:raise ValueError('Failed to bind rank')
    os.environ['WAVETRISK_EXPECT_CPU200']=str(cpu)
    os.environ['MV2_ENABLE_AFFINITY']='0'
    if sys.argv[2]=='--probe':
        print(json.dumps(placement(rank,[])),flush=True);return 0
    return execute(Path(sys.argv[2]),['numactl','--localalloc',*sys.argv[3:]])
if __name__=='__main__':raise SystemExit(main())
