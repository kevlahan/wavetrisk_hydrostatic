#!/usr/bin/env python3
"""Pin one rank to its planned physical core before starting MPI."""
import json
import os
from pathlib import Path
import socket
import sys
from rank_launch_190 import execute, placement


def pin(rank, cpu):
    allowed = sorted(os.sched_getaffinity(0))
    if cpu not in allowed:
        raise ValueError(f'Rank {rank} on {socket.gethostname()}: requested CPU {cpu} '
                         f'unavailable; allowed CPUs={allowed}. Check Slurm step allocation.')
    os.sched_setaffinity(0, {cpu})
    os.environ['WAVETRISK_EXPECT_CPU200'] = str(cpu)
    os.environ['MV2_ENABLE_AFFINITY'] = '0'


def main():
    case = json.loads(Path(sys.argv[1]).read_text())
    rank = int(os.environ['SLURM_PROCID'])
    pin(rank, case['mapping'][rank]['cpu'])
    if sys.argv[2] == '--probe':
        print(json.dumps(placement(rank, [])), flush=True)
        return 0
    return execute(Path(sys.argv[2]), ['numactl', '--localalloc', *sys.argv[3:]])


if __name__ == '__main__':
    raise SystemExit(main())
