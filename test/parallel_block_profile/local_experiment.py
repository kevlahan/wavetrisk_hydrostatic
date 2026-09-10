#!/usr/bin/env python3
"""Sequential four-rank local diagnosis, NOT an 83-rank performance forecast.

Frozen fresh inputs, alternating order, no build subprocesses. Source-instrumented
legacy remains explicitly separate from authoritative legacy. No core-placement
claim on macOS. Optional sampling belongs in a separate diagnostic execution.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import time

from analyze import analyze
from experiment import SWITCHES, digest


def legacy_detail(directory, ranks):
    rows = []
    for rank in range(ranks):
        steps = {}
        counts, boundary = {}, {}
        for line in (directory / f'legacy-detail-rank-{rank}.txt').read_text().splitlines():
            f = line.split()
            if not f:
                continue
            if f[0] == 'step':
                steps[int(f[1])] = {'restart': int(f[2]), 'remap': int(f[3]),
                                    'wall': float(f[4]), 'cpu': float(f[5]), 'self': {}}
            elif f[0] == 'self':
                steps[int(f[1])]['self'][int(f[2])] = [float(f[3]), float(f[4])]
            elif f[0] == 'counter':
                counts[int(f[1])] = int(f[2])
            elif f[0] == 'boundary':
                boundary[int(f[1])] = [float(f[2]), float(f[3]), int(f[4])]
        for step in steps.values():
            if len(step['self']) != 65:
                raise ValueError('Incomplete legacy profile')
            if abs(sum(v[0] for v in step['self'].values()) - step['wall']) > 1e-6:
                raise ValueError('Legacy self wall conservation failed')
            if abs(sum(v[1] for v in step['self'].values()) - step['cpu']) > 1e-6:
                raise ValueError('Legacy self CPU conservation failed')
        rows.append({'rank': rank, 'steps': steps, 'counters': counts, 'boundary': boundary})
    return rows


def run(args):
    if args.pairs:
        raise ValueError('Printed-state-only timing is retired. Use local_campaign.py for mandatory full-field validation.')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    seed = args.fixture.resolve(strict=True)
    binaries = {'legacy': args.legacy.resolve(strict=True), 'block': args.block.resolve(strict=True),
                'legacy-off': args.instrumented.resolve(strict=True),
                'legacy-detail': args.instrumented.resolve(strict=True), 'block-detail': args.block.resolve(strict=True)}
    hashes = {k: digest(v) for k, v in binaries.items()}
    input_text = (seed / 'simple.in').read_text()
    input_text, count = re.subn(r'(?m)^time_end\s+\S+', f'time_end           {args.time_end}', input_text)
    if count != 1:
        raise ValueError('Expected one time_end')
    checkpoint = seed / 'test_checkpoint_0003.bin.zst'
    checkpoint_hash = digest(checkpoint)
    schedule = []
    for pair in range(args.pairs):
        order = ['legacy', 'block'] if pair % 2 == 0 else ['block', 'legacy']
        schedule.extend((pair, label) for label in order)
    if args.detail:
        schedule += [(0, 'legacy-off'), (0, 'legacy-detail'), (0, 'block-detail')]
    if args.only:
        if args.pairs:
            raise ValueError('--only requires --pairs 0 (diagnosis, not paired timing)')
        labels = args.only.split(',')
        if any(label not in binaries for label in labels):
            raise ValueError('Unknown --only label')
        schedule = list(enumerate(labels))
    records = []
    reference = None
    env = {k: v for k, v in os.environ.items() if not k.startswith('WAVETRISK_')}
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
        env[key] = '1'
    for index, (pair, label) in enumerate(schedule):
        path = out / f'{index:02d}-{label}'
        path.mkdir()
        (path / 'simple.in').write_text(input_text)
        shutil.copy2(checkpoint, path / checkpoint.name)
        (path / 'grids').symlink_to((seed / 'grids').resolve(), target_is_directory=True)
        shutil.copy2(binaries[label], path / 'climateJ5')
        if digest(path / 'climateJ5') != hashes[label] or digest(checkpoint) != checkpoint_hash:
            raise ValueError('Binary/input changed')
        for key in SWITCHES:
            env[key] = '0'
        if label.endswith('-detail'):
            env['WAVETRISK_PROFILE_BLOCK_DETAIL'] = '1'
            env['WAVETRISK_PROFILE_PARALLEL_BLOCKS'] = '1'
        command = ['mpirun', '-n', str(args.ranks), './climateJ5', 'simple.in']
        print(f'Start {index + 1}/{len(schedule)} {label}: {path}', flush=True)
        start = time.monotonic()
        with (path / 'run.log').open('x') as log:
            result = subprocess.run(command, cwd=path, env=env, stdout=log, stderr=subprocess.STDOUT)
        elapsed = time.monotonic() - start
        if result.returncode:
            raise RuntimeError(f'{label} exit {result.returncode}; see {path / "run.log"}')
        report = analyze(path / 'run.log')
        if not report['non_checkpoint_step_times']:
            raise ValueError('No measured steps')
        # Exact agreement at printed precision, NOT a full-field oracle.
        if not report['completed']:
            raise ValueError(f'Incomplete/failed run: {path}')
        state = report['printed_states']
        if reference is None:
            reference = state
        if state != reference:
            raise ValueError(f'Printed states differ: {path}')
        record = {'label': label, 'pair': pair, 'directory': str(path), 'launch_elapsed': elapsed,
                  'binary_sha256': hashes[label], 'analysis': report, 'environment': {k: env[k] for k in SWITCHES}}
        if label == 'legacy-detail':
            record['legacy_detail'] = legacy_detail(path, args.ranks)
        cp = path / 'test_checkpoint_0004.bin.zst'
        if cp.exists():
            record['checkpoint4_sha256'] = digest(cp)
        records.append(record)
        output = {'host': platform.uname()._asdict(), 'ranks': args.ranks, 'placement_verified': False,
                  'fixture_checkpoint_sha256': checkpoint_hash, 'records': records}
        (out / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
        print(f'Completed {label}: steps {report["non_checkpoint_step_times"]}', flush=True)
    ratios = []
    for pair in range(args.pairs):
        times = {r['label']: statistics.mean(r['analysis']['non_checkpoint_step_times'])
                 for r in records if r['pair'] == pair and r['label'] in ('legacy', 'block')}
        ratios.append(times['block'] / times['legacy'])
    print('Paired mean-step ratios:', ratios, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('legacy', 'instrumented', 'block', 'fixture', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--pairs', type=int, default=3)
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--time-end', default='0.3040')
    parser.add_argument('--detail', action='store_true')
    parser.add_argument('--only', help='Comma-separated labels for a separate diagnostic run; requires --pairs 0')
    run(parser.parse_args())
