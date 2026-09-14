#!/usr/bin/env python3
"""Record rank placement and child resources using wait4, without GNU time.

Each rank's Python parent sleeps while the solver runs. RSS/fault counts describe
that solver child, not the parent or srun. No sampling or solver instrumentation.
"""
import json
import os
from pathlib import Path
import socket
import subprocess
import sys


def placement(rank, command):
    status=Path('/proc/self/status')
    cpuinfo=Path('/proc/cpuinfo')
    return {'rank':rank,'hostname':socket.gethostname(),'launcher_pid':os.getpid(),
            'affinity':sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
            'command':command,'resource_collector':'os.wait4',
            'resource_scope':'Solver child peak RSS and faults; excludes the waiting Python parent and MPI launcher.',
            'memory_affinity':next((line for line in status.read_text().splitlines()
                                   if line.startswith('Mems_allowed_list:')),'unavailable') if status.exists() else 'unavailable',
            'cpu_model':next((line.split(':',1)[1].strip() for line in cpuinfo.read_text().splitlines()
                              if line.startswith('model name')),'unavailable') if cpuinfo.exists() else 'unavailable',
            'slurm':{k:os.environ[k] for k in ('SLURM_JOB_ID','SLURM_LOCALID','SLURM_NODEID') if k in os.environ}}


def execute(directory,command):
    if not command or 'SLURM_PROCID' not in os.environ:
        raise ValueError('Use under srun: rank_launch_190.py DIRECTORY executable [arguments]')
    if not hasattr(os,'wait4'): raise RuntimeError('This runner requires Unix os.wait4')
    rank=int(os.environ['SLURM_PROCID'])
    result=placement(rank,command)
    path=directory/f'rank-{rank}.json'
    with path.open('x') as stream: json.dump(result,stream,indent=2)
    try:
        # Preserve descriptors supplied by MPI/PMI launchers, as exec would.
        process=subprocess.Popen(command,close_fds=False)
    except OSError as error:
        result.update(returncode=127,launch_error=str(error))
        path.write_text(json.dumps(result,indent=2)+'\n')
        print(f'Rank {rank}: cannot launch {command[0]}: {error}',file=sys.stderr,flush=True)
        return 127
    while True:
        try:
            pid,status,usage=os.wait4(process.pid,0)
            break
        except InterruptedError:
            continue
    process.returncode=os.waitstatus_to_exitcode(status)
    # Linux reports KiB; Darwin reports bytes. The cluster protocol uses Linux.
    rss=usage.ru_maxrss/1024 if sys.platform=='darwin' else usage.ru_maxrss
    result.update(solver_pid=pid,returncode=process.returncode,max_rss_kib=int(rss),
                  major_faults=usage.ru_majflt,minor_faults=usage.ru_minflt,
                  user_cpu_seconds=usage.ru_utime,system_cpu_seconds=usage.ru_stime)
    path.write_text(json.dumps(result,indent=2)+'\n')
    # Retain the three existing parser labels, explicitly identifying the collector.
    path.with_suffix('.resources').write_text(
        'Resource collector: Python os.wait4 (solver child)\n'
        f'Maximum resident set size (kbytes): {int(rss)}\n'
        f'Major (requiring I/O) page faults: {usage.ru_majflt}\n'
        f'Minor (reclaiming a frame) page faults: {usage.ru_minflt}\n')
    return process.returncode if process.returncode>=0 else 128-process.returncode


if __name__=='__main__':
    raise SystemExit(execute(Path(sys.argv[1]),sys.argv[2:]))
