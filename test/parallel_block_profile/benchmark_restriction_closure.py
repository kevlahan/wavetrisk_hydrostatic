#!/usr/bin/env python3
"""Time actual scalar and column restriction closure on a bounded synthetic storage fixture.

Run sequentially on an idle machine, outside solver timing and builds. This
isolates restriction closure/storage access; it does not model MPI or whole-run memory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import threading
import time
from test_restriction_closure import ROOT, prepare
from profile_protocol import MemoryMonitor


from benchmark_restriction_columns import ComponentMonitor


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    prepare(out)
    command = ['gfortran','-cpp', '-O2', '-march=native', '-ffp-contract=off', '-I' + str(out),
               str(ROOT / 'src/kind.f90'), str(ROOT / 'src/parallel_block_scalar_storage.f90'),
               'closure.f90', '-o', 'benchmark']
    built = subprocess.run(command, cwd=out, capture_output=True, text=True)
    (out / 'build.log').write_text(built.stdout + built.stderr)
    built.check_returncode()
    monitor = ComponentMonitor(out)
    with (out / 'run.log').open('w') as log:
        measured = subprocess.Popen(['./benchmark', 'benchmark'], cwd=out, stdout=log, stderr=subprocess.STDOUT)
        monitor.start(measured.pid)
        code = measured.wait()
    memory = monitor.finish()
    memory['scope'] = 'System-wide VM deltas and sampled benchmark-process RSS; excludes compressed footprint'
    (out / 'memory.json').write_text(json.dumps(memory, indent=2) + '\n')
    if code:
        raise subprocess.CalledProcessError(code, './benchmark benchmark')
    rows = []
    for line in (out / 'run.log').read_text().splitlines():
        fields = line.split()
        if len(fields) == 5:
            rows.append(dict(trial=int(fields[0]), variant=int(fields[1]),
                             cpu_seconds=float(fields[2]), wall_seconds=float(fields[3]),
                             checksum=float(fields[4])))
    if len(rows) != 8 or len({r['checksum'] for r in rows}) != 1:
        raise ValueError('Incomplete benchmark or mismatched checksums')
    medians = {name: statistics.median(r['wall_seconds'] for r in rows if r['variant'] == variant)
               for variant, name in [(0, 'scalar'), (1, 'column')]}
    result = {'build_command': command,
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [ROOT / 'src/kind.f90', ROOT / 'src/parallel_block_scalar_storage.f90',
                                          *sorted(out.glob('*.inc')), out / 'closure.f90']},
              'memory': memory, 'trials': rows, 'median_wall_seconds': medians,
              'reduction_percent': 100 * (1 - medians['column'] / medians['scalar']),
              'scope': 'Synthetic component fixture; no whole-solver or cluster speedup claim.'}
    (out / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    run(parser.parse_args().out.resolve())
