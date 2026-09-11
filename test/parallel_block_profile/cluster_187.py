#!/usr/bin/env python3
"""Sequential J5/J7 cluster restart gates and baseline/candidate timing.

Run inside an 83-task Slurm allocation. Requires Python 3.9+, numpy, zstd and
srun. No compilation, submission, source modification or tolerance adjustment.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time

from analyze import analyze, STEP
from checkpoint_gate import accept
from compare_restart_checkpoints import compare
from profile_protocol import parameters, checkpoint_name, grid_identity, environment


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def copy_grids(source, destination, expected=None):
    required = source / 'HR95JT_WT_004'
    if not required.is_file() or required.stat().st_size == 0:
        raise ValueError('Missing/empty J5 optimized grid file: ' + str(required))
    # Materialize directory and file symlinks before any measured MPI work.
    shutil.copytree(source, destination, symlinks=False)
    actual = grid_identity(destination)
    if expected is not None and actual != expected:
        raise ValueError('Grid copy differs from the frozen input snapshot')
    return actual


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def schedule(repeats, profiles):
    jobs = [('checked-off', 'checked', False, False, 'validation'),
            ('checked-on', 'checked', True, False, 'validation')]
    jobs += [(kind + '-warmup', kind, False, False, 'warmup')
             for kind in ('baseline', 'candidate')]
    for i in range(1, repeats + 1):
        # Reverse each pair's order to reduce monotonic placement/load drift.
        kinds = ('baseline', 'candidate') if i % 2 else ('candidate', 'baseline')
        jobs += [(f'{kind}-{i}', kind, False, False, 'timing') for kind in kinds]
    if profiles:
        jobs += [(kind + '-profile', kind, False, True, 'profile')
                 for kind in ('baseline', 'candidate')]
    return jobs


def validate_log(root, run_id, seed_index, expected, require_remap=False):
    log = (root / 'run.log').read_text(errors='replace')
    a = analyze(root / 'run.log')
    if not a['completed'] or not a['step_records']:
        raise ValueError('Incomplete run or runtime failure: ' + str(root))
    for pattern in (r'min_level\s*=\s*5\b', r'DOMAIN_LEVEL\s*=\s*2\b'):
        if not re.search(pattern, log):
            raise ValueError('Log does not identify a param_J5 build: ' + pattern)
    checkpoints = {p.name for p in root.glob('*_checkpoint_*.bin.zst')}
    wanted = {f'{run_id}_checkpoint_{i:04d}.bin.zst' for i in (seed_index, *expected)}
    if checkpoints != wanted:
        raise ValueError('Unexpected checkpoint set: ' + repr(sorted(checkpoints)))
    reloads = list(map(int, re.findall(r'Restarting from checkpoint\s+(\d+)', log)))
    if reloads != [seed_index, *expected]:
        raise ValueError('Missing or unexpected reloads: ' + repr(reloads))
    post_restart_remaps = {}
    for index in expected:
        match = re.search(r'Restarting from checkpoint\s+' + str(index) + r'\b', log)
        after = log[match.end():]
        if not any(STEP.match(line) for line in after.splitlines()):
            raise ValueError('No timestep after reload ' + str(index))
        post_restart_remaps[str(index)] = 'Remapping vertical coordinates' in after
        if require_remap and not post_restart_remaps[str(index)]:
            raise ValueError('No remap after reload ' + str(index))
    return a, {'reloads': reloads, 'post_restart_remaps': post_restart_remaps}


def exact_gate(left, right, run_id, expected):
    reports, failures = {}, []
    for index in expected:
        name = f'{run_id}_checkpoint_{index:04d}.bin.zst'
        report = compare(left / name, right / name, expected_domains=160)
        reports[name] = report
        failures.extend(name + ': ' + reason for reason in accept(report, exact=True))
    return {'passed': not failures, 'failures': failures, 'checkpoints': reports,
            'scope': 'Exact semantic fields, thresholds and topology; compressed bytes/load weights may differ.'}


def summarize(runs):
    result = {}
    for kind in ('baseline', 'candidate'):
        selected = [r for r in runs.values() if r['kind'] == kind and r['phase'] == 'timing']
        elapsed = [r['elapsed_seconds'] for r in selected]
        steps = [r['ordinary_step_median'] for r in selected]
        result[kind] = {'elapsed_seconds': elapsed, 'ordinary_step_medians': steps,
                        'median_elapsed_seconds': statistics.median(elapsed),
                        'median_ordinary_step_seconds': statistics.median(steps)}
    result['observed_elapsed_reduction_percent'] = 100 * (
        1 - result['candidate']['median_elapsed_seconds'] / result['baseline']['median_elapsed_seconds'])
    result['scope'] = ('Observed timings only. Inspect node load, placement and memory/swap evidence. '
                       'Warmups, correctness runs and detailed profiles are excluded.')
    return result


def run(args):
    fixture = args.fixture.resolve(strict=True)
    out = args.out.resolve()
    binaries = {kind: getattr(args, kind).resolve(strict=True) for kind in ('baseline', 'candidate', 'checked')}
    if out.exists() or out.is_relative_to(fixture) or fixture.is_relative_to(out):
        raise ValueError('Output must be a new sibling directory outside the fixture')
    for binary in binaries.values():
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError('Not an executable: ' + str(binary))
        if binary.is_relative_to(out):
            raise ValueError('Output overlaps an executable')
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Run inside a Slurm allocation (salloc or sbatch)')
    allocation_tasks = os.environ.get('SLURM_NTASKS')
    if allocation_tasks and int(allocation_tasks) < args.ranks:
        raise ValueError('Allocation has fewer tasks than requested')
    for command in ('srun', 'zstd'):
        if not shutil.which(command):
            raise ValueError('Missing command: ' + command)
    if args.ranks != 83 or args.repeats < 1:
        raise ValueError('This protocol requires 83 ranks and at least one measured pair')
    text = (fixture / 'simple.in').read_text()
    p = parameters(text)
    for key, value in {'max_level': '7', 'zlevels': '30', 'Nsoil': '10',
                       'physics_type': 'Simple', 'NCAR_topo': 'F', 'sso': 'F',
                       'resume': '3', 'CP_EVERY': '1'}.items():
        if p.get(key) != value:
            raise ValueError(f'Expected {key}={value}; got {p.get(key)}')
    if float(p['time_end']) <= 0.3 or float(p['dt_write']) <= 0:
        raise ValueError('Invalid restart endpoint/write interval')
    expected = args.expected_checkpoints
    if expected != sorted(set(expected)) or not expected or expected[0] <= 3:
        raise ValueError('Expected new checkpoints must be unique, increasing and greater than 3')
    seed = checkpoint_name(text)
    # Cheap preflight verifies the header without decompressing the large seed.
    import struct
    with subprocess.Popen(['zstd', '-dc', str(fixture / seed)], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL) as process:
        header = process.stdout.read(24)
        process.stdout.close()
        process.wait()
    if len(header) != 24 or struct.unpack('<3q', header) != (0x5741564554524953, 1, 160):
        raise ValueError('Seed must be a genuine J5 checkpoint with 160 Domains')
    if not (fixture / 'grids' / 'HR95JT_WT_004').is_file():
        raise ValueError('Fixture lacks grids/HR95JT_WT_004: ' + str(fixture))
    out.mkdir(parents=True)
    inputs = out / 'inputs'
    inputs.mkdir()
    grids = copy_grids(fixture / 'grids', inputs / 'grids')
    shutil.copy2(fixture / 'simple.in', inputs / 'simple.in')
    shutil.copy2(fixture / seed, inputs / seed)
    for kind, binary in binaries.items():
        shutil.copy2(binary, inputs / kind)
    jobs = schedule(args.repeats, args.profiles)
    result = {'status': 'running', 'parameters': p, 'ranks': args.ranks,
              'fixture': str(fixture), 'input_sha256': digest(inputs / 'simple.in'),
              'seed_sha256': digest(inputs / seed), 'grids': grids,
              'executables': {k: {'path': str(v), 'sha256': digest(inputs / k)} for k, v in binaries.items()},
              'build_notes': args.build_notes, 'slurm_environment': {
                  k: v for k, v in os.environ.items() if k.startswith('SLURM_')},
              'schedule': [j[0] for j in jobs], 'runs': {}, 'gates': {},
              'limitations': ['No checked-versus-optimized or legacy numerical equivalence claim.',
                             'No automatic distributed memory-pressure qualification.',
                             'Checked solver flags do not necessarily cover external physics.']}
    result_path = out / 'results.json'
    save(result_path, result)
    try:
        for label, kind, oracle, profile, phase in jobs:
            result['active'] = label
            save(result_path, result)
            print('Starting ' + label, flush=True)
            directory = out / label
            directory.mkdir()
            for name in ('simple.in', seed):
                shutil.copy2(inputs / name, directory / name)
            shutil.copy2(inputs / kind, directory / 'climateJ5')
            copy_grids(inputs / 'grids', directory / 'grids', expected=grids)
            env = environment(oracles=oracle, detail=profile)
            command = ['srun', '--ntasks=83', '--cpus-per-task=1', '--cpu-bind=cores',
                       '--kill-on-bad-exit=1', './climateJ5', 'simple.in']
            start = time.monotonic()
            with (directory / 'run.log').open('x') as log:
                code = subprocess.run(command, cwd=directory, env=env,
                                      stdout=log, stderr=subprocess.STDOUT).returncode
            elapsed = time.monotonic() - start
            record = {'kind': kind, 'phase': phase, 'returncode': code, 'elapsed_seconds': elapsed,
                      'command': command, 'oracles': oracle, 'profile': profile,
                      'binary_sha256': digest(directory / 'climateJ5'),
                      'environment': {k: v for k, v in env.items() if k.startswith('WAVETRISK_') or
                                      k.endswith('NUM_THREADS') or k == 'VECLIB_MAXIMUM_THREADS'}}
            result['runs'][label] = record
            save(result_path, result)
            if code:
                raise RuntimeError(f'{label} exited {code}; see {directory / "run.log"}')
            analysis, coverage = validate_log(directory, p['run_id'], 3, expected, args.require_post_restart_remap)
            record.update(coverage, step_count=len(analysis['step_records']),
                          ordinary_step_median=analysis['non_checkpoint_median'])
            save(directory / 'analysis.json', analysis)
            if record['ordinary_step_median'] is None:
                raise ValueError('No ordinary steps available for comparison')
            # Self-comparisons of the first references check finiteness/layout only.
            # All later comparisons use an independent run of the matching build mode.
            reference = out / ('checked-off' if phase == 'validation' else 'baseline-warmup')
            gate = exact_gate(reference, directory, p['run_id'], expected)
            gate['reference'] = str(reference)
            gate['self_check'] = reference == directory
            result['gates'][label] = gate
            save(result_path, result)
            if not gate['passed']:
                raise RuntimeError(f'Exact checkpoint gate failed for {label}: {gate["failures"]}')
            print(f'PASS {label}: {elapsed:.2f}s; ordinary-step median {record["ordinary_step_median"]:.4g}s', flush=True)
        if (grid_identity(inputs / 'grids') != grids or digest(inputs / seed) != result['seed_sha256'] or
                digest(inputs / 'simple.in') != result['input_sha256']):
            raise ValueError('Frozen input snapshot changed during the campaign')
        result['summary'] = summarize(result['runs'])
        result['status'] = 'passed'
        result.pop('active', None)
        save(result_path, result)
        print(json.dumps(result['summary'], indent=2), flush=True)
    except BaseException as error:
        result['status'] = 'failed'
        result['error'] = repr(error)
        save(result_path, result)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'candidate', 'checked', 'fixture', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=83)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--expected-checkpoints', type=int, nargs='+', default=[4])
    parser.add_argument('--require-post-restart-remap', action='store_true')
    parser.add_argument('--profiles', action='store_true', help='Append separate baseline and candidate profiles')
    parser.add_argument('--build-notes', required=True, help='Compiler/MPI versions, flags and source identifiers')
    run(parser.parse_args())


if __name__ == '__main__':
    main()
