#!/usr/bin/env python3
"""Separate macOS call-stack diagnostic; NEVER include these runs in timing pairs.

Samples rank 0 only, beginning at the application simulation marker. macOS sample
includes sleeping stacks; sampled percentages are not useful-CPU percentages.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from experiment import SWITCHES, digest

if len(sys.argv) > 1 and sys.argv[1] == '--rank-exec':
    rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
    Path(f'rank-{rank}.json').write_text(json.dumps({'rank': rank, 'pid': os.getpid()}))
    os.execv('./climateJ5', ['./climateJ5', 'simple.in'])

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--legacy', type=Path, required=True)
parser.add_argument('--block', type=Path, required=True)
parser.add_argument('--fixture', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)
env = {k: v for k, v in os.environ.items() if not k.startswith('WAVETRISK_')}
for key in SWITCHES:
    env[key] = '0'
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    env[key] = '1'
for label, binary in [('legacy', args.legacy), ('block', args.block)]:
    directory = args.out.resolve() / label
    directory.mkdir()
    shutil.copy2(binary, directory / 'climateJ5')
    shutil.copy2(args.fixture / 'test_checkpoint_0003.bin.zst', directory)
    (directory / 'grids').symlink_to((args.fixture / 'grids').resolve(), target_is_directory=True)
    text = re.sub(r'(?m)^time_end\s+\S+', 'time_end           0.3040', (args.fixture / 'simple.in').read_text())
    (directory / 'simple.in').write_text(text)
    command = ['mpirun', '-n', '4', sys.executable, str(Path(__file__).resolve()), '--rank-exec']
    print(f'Start separate sampled {label}: {directory}', flush=True)
    with (directory / 'run.log').open('x') as log:
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT)
        sampler = None
        while process.poll() is None:
            if sampler is None and 'Start simulation run' in (directory / 'run.log').read_text():
                identity = json.loads((directory / 'rank-0.json').read_text())
                sampler = subprocess.Popen(['/usr/bin/sample', str(identity['pid']), '20', '10',
                                            '-mayDie', '-fullPaths', '-file', str(directory / 'rank-0.sample.txt')],
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                print(f'Sampling {label} rank 0 pid {identity["pid"]}', flush=True)
            time.sleep(1)
        if sampler:
            output, _ = sampler.communicate()
            (directory / 'sample-status.txt').write_text(output)
        if process.returncode:
            raise RuntimeError(f'Failed {label}: {process.returncode}')
    (directory / 'identity.json').write_text(json.dumps({'binary_sha256': digest(binary),
             'diagnostic_only': True, 'command': command, 'sampler_exit': sampler.returncode if sampler else None}, indent=2))
    print(f'Completed sampled {label}', flush=True)
