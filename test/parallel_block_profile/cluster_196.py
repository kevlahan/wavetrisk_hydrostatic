#!/usr/bin/env python3
"""Stage 196: paired legacy/block workflow diagnosis with observer calibration.

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
import sys
import subprocess
import time

from analyze import analyze, STEP
from legacy_metrics_191 import collect as collect_legacy
from rank_metrics_190 import snapshot_allocation, node_memory, rank_resources
from checkpoint_gate import accept, POLICY
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
    jobs=[('checked-off','checked',False,False,'validation'),('checked-on','checked',True,False,'validation'),
          ('baseline-warmup','baseline',False,False,'warmup'),('candidate-warmup','candidate',False,False,'warmup'),
          ('legacy-warmup','legacy',False,False,'warmup'),('legacy-minimal-warmup','legacy_profile',False,False,'warmup')]
    for mode in ('coarse','detail'):
        for side,kind in [('block','candidate'),('legacy','legacy_profile')]:
            jobs.append((side+'-'+mode+'-warmup',kind,False,False,'warmup'))
    modes=('bare','minimal','coarse','detail')
    for i in range(repeats):
        for mode in modes[i%4:]+modes[:i%4]:
            sides=('block','legacy') if i%2==0 else ('legacy','block')
            for side in sides:
                kind=('baseline' if side=='block' else 'legacy') if mode=='bare' else ('candidate' if side=='block' else 'legacy_profile')
                jobs.append((f'{side}-{mode}-{i+1}',kind,False,False,'timing-'+mode))
    return jobs

def observer_mode(label,phase):
    return 1

def operation_mode(label,phase):
    if phase=='validation' or '-detail-' in label:return 2
    if '-coarse-' in label:return 1
    return 0


def stage191_evidence(directory, analysis, install, overlap, profile):
    text=(directory/'run.log').read_text()
    matches=re.findall(r'Stage 191 requested install/mass-overlap =\s*([01])\s+([01])',text)
    if matches != [(install,overlap)]:
        raise ValueError('Missing or inconsistent Stage 191 selector record')
    if profile:
        windows=analysis['detail_windows']
        required=('source proof payload bytes','receive proof payload bytes','copy proof payload bytes',
                  'install copied payload bytes','install moved payload bytes')
        totals={key:0 for key in required}
        if not windows: raise ValueError('Missing Stage 191 installation profile')
        for window in windows:
            counters={k.strip():v for k,v in window['counters'].items()}
            if any(k not in counters for k in required): raise ValueError('Missing installation counters')
            for key in totals: totals[key]+=counters[key]
        if any(totals[key] for key in required[:-1]) or totals[required[-1]]<=0:
            raise ValueError('Production installation did not move storage without proofs/copies')
        return totals
    return {'requested_install':install,'requested_overlap':overlap}


def stage194_evidence(directory, columns):
    matches=re.findall(r'Stage 194 requested inner inverse geometry = ([01])', (directory/'run.log').read_text())
    if matches != [columns]: raise ValueError('Missing or inconsistent Stage 194 selector record')
    return {'requested_inner_geometry':columns}


def legacy_gate(left, right, run_id, expected):
    reports, failures = {}, []
    for index in expected:
        name = f'{run_id}_checkpoint_{index:04d}.bin.zst'
        report = compare(left/name, right/name, expected_domains=160)
        reports[name] = report
        failures.extend(name+': '+reason for reason in accept(report,exact=False))
    return {'passed':not failures,'failures':failures,'checkpoints':reports,'policy':POLICY,
            'scope':'Existing fixed numerical screen; not bitwise legacy equivalence.'}


def workload_gate(reference, candidate):
    # State records exclude CPU timings and include step/time/dt/Jmax/dof.
    left=[r['state'] for r in reference['step_records']]
    right=[r['state'] for r in candidate['step_records']]
    return {'passed':left==right,'reference_steps':left,'candidate_steps':right,
            'scope':'Identical printed timestep, dt, Jmax and dof sequence required for timing.'}


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


from analyze_workflow_195 import collect as collect_workflow
from analyze_operations_196 import collect as collect_operations, summarize


def read_build_notes(args):
    path = getattr(args, 'build_notes_file', None)
    if path is None:
        return args.build_notes
    # Preserve the complete provenance in results.json without putting it in argv.
    notes = path.read_text()
    if not isinstance(json.loads(notes), dict):
        raise ValueError('Build notes file must contain a JSON object')
    return notes


def run(args):
    build_notes = read_build_notes(args)
    fixture = args.fixture.resolve(strict=True)
    out = args.out.resolve()
    binaries = {kind: getattr(args, kind).resolve(strict=True) for kind in ('baseline', 'candidate', 'checked', 'legacy', 'legacy_profile')}
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
    placement = snapshot_allocation(out, args.ranks)
    result = {'status': 'running', 'parameters': p, 'ranks': args.ranks,
              'fixture': str(fixture), 'input_sha256': digest(inputs / 'simple.in'),
              'seed_sha256': digest(inputs / seed), 'grids': grids,
              'executables': {k: {'path': str(v), 'sha256': digest(inputs / k)} for k, v in binaries.items()},
              'build_notes': build_notes, 'slurm_environment': {
                  k: v for k, v in os.environ.items() if k.startswith('SLURM_')},
              'placement': placement, 'schedule': [j[0] for j in jobs], 'runs': {}, 'gates': {},
              'limitations': ['Legacy Domain numerics plus common I/O fix; fixed numerical screen, not bitwise legacy equality.',
                             'Per-rank RSS and NUMA placement snapshots are diagnostic; node VM counters cover the whole allocation.',
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
            install,overlap = ('1','0')
            columns = '1'
            env.update(WAVETRISK_STAGE191_INSTALL=install,WAVETRISK_STAGE191_MASS_OVERLAP=overlap,
                       WAVETRISK_STAGE192_SPARSE_OUTER='0',WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY=columns)
            env['WAVETRISK_WORKFLOW195']='1'
            env['WAVETRISK_OPERATIONS196']=str(operation_mode(label,phase))
            command = ['srun', '--ntasks=83', '--cpus-per-task=1', '--cpu-bind=cores',
                       '--kill-on-bad-exit=1', sys.executable, str(Path(__file__).resolve().with_name('rank_launch_190.py')), str(directory), './climateJ5', 'simple.in']
            memory_before = node_memory(out / (label+'-memory-before.json'))
            start = time.monotonic()
            with (directory / 'run.log').open('x') as log:
                code = subprocess.run(command, cwd=directory, env=env,
                                      stdout=log, stderr=subprocess.STDOUT).returncode
            elapsed = time.monotonic() - start
            memory_after = node_memory(out / (label+'-memory-after.json'))
            record = {'side':'legacy' if kind in ('legacy','legacy_profile') else 'block',
                      'round':int(label.rsplit('-',1)[1]) if phase.startswith('timing-') else None, 'kind': kind, 'phase': phase, 'returncode': code, 'elapsed_seconds': elapsed,
                      'command': command, 'oracles': oracle, 'profile': profile,
                      'binary_sha256': digest(directory / 'climateJ5'),
                      'memory_before': memory_before, 'memory_after': memory_after,
                      'rank_metrics': [],
                      'environment': {k: v for k, v in env.items() if k.startswith('WAVETRISK_') or
                                      k.endswith('NUM_THREADS') or k == 'VECLIB_MAXIMUM_THREADS'}}
            result['runs'][label] = record
            save(result_path, result)
            if code:
                raise RuntimeError(f'{label} exited {code}; see {directory / "run.log"}')
            # A failed launcher may never create child resource files.
            # Preserve its exit status and log instead of masking the cause.
            record['rank_metrics'] = rank_resources(directory)
            if len(record['rank_metrics']) != args.ranks:
                raise ValueError('Missing per-rank exit/RSS records: '+label)
            save(result_path, result)
            analysis, coverage = validate_log(directory, p['run_id'], 3, expected, args.require_post_restart_remap)
            if kind in ('baseline','candidate','checked'):
                if re.findall(r'Stage 193 requested scalar inverse columns = ([01])',(directory/'run.log').read_text())!=['1']:
                    raise ValueError('Stage 193 scalar kernel must be enabled')
            if kind in ('candidate','checked'):
                record['stage191_evidence']=stage191_evidence(directory,analysis,install,overlap,profile)
                record['stage194_evidence']=stage194_evidence(directory,columns)
            record.update(coverage, step_count=len(analysis['step_records']),
                          ordinary_step_median=analysis['non_checkpoint_median'])
            save(directory / 'analysis.json', analysis)
            if record['ordinary_step_median'] is None:
                raise ValueError('No ordinary steps available for comparison')
            # Self-comparisons of the first references check finiteness/layout only.
            # All later comparisons use an independent run of the matching build mode.
            reference = out / ('checked-off' if phase == 'validation' else
                               'legacy-warmup' if kind in ('legacy','legacy_profile') else 'baseline-warmup')
            gate = exact_gate(reference, directory, p['run_id'], expected)
            gate['reference'] = str(reference)
            gate['self_check'] = reference == directory
            result['gates'][label] = gate
            save(result_path, result)
            if not gate['passed']:
                raise RuntimeError(f'Exact checkpoint gate failed for {label}: {gate["failures"]}')
            if phase != 'validation':
                reference_analysis = json.loads((out/'baseline-warmup'/'analysis.json').read_text())
                work = workload_gate(reference_analysis, analysis)
                result['gates']['workload-'+label] = work
                if not work['passed']:
                    raise ValueError('Workload differs; timing comparison rejected: '+label)
            if label=='baseline-warmup':
                for checked_label in ('checked-off','checked-on'):
                    checked_analysis=json.loads((out/checked_label/'analysis.json').read_text())
                    work=workload_gate(checked_analysis,analysis)
                    result['gates']['workload-'+checked_label]=work
                    if not work['passed']:raise ValueError('Checked workload differs from reference')
            if label in ('baseline-warmup','candidate-warmup'):
                gate = exact_gate(out/'checked-off',directory,p['run_id'],expected)
                result['gates']['cross-mode-'+label] = gate
                if not gate['passed']:
                    raise ValueError('Checked/optimized exact mismatch: '+label)
            if kind in ('legacy','legacy_profile'):
                gate = legacy_gate(directory,out/'candidate-warmup',p['run_id'],expected)
                result['gates']['legacy-screen-'+label] = gate
                if not gate['passed']:
                    raise ValueError('Legacy numerical screen failed: '+label)
            if kind=='legacy_profile' and profile: record['legacy_profile']=collect_legacy(directory,args.ranks)
            if kind in ('candidate','checked','legacy_profile'):
                record['workflow']=collect_workflow(directory,args.ranks,observer_mode(label,phase))
                record['operations']=collect_operations(directory,args.ranks,operation_mode(label,phase),record['side'])
                if len(record['workflow']['steps'])!=record['step_count']:raise ValueError('Workflow coverage differs')
                if sum(s['restart'] for s in record['workflow']['steps'])!=len(expected):raise ValueError('Workflow restart flags differ')
                if sum(s['remap'] for s in record['workflow']['steps'])!=(directory/'run.log').read_text().count('Remapping vertical coordinates'):raise ValueError('Workflow remap flags differ')
            save(result_path,result)
            print(f'PASS {label}: {elapsed:.2f}s; ordinary-step median {record["ordinary_step_median"]:.4g}s', flush=True)
        if (grid_identity(inputs / 'grids') != grids or digest(inputs / seed) != result['seed_sha256'] or
                digest(inputs / 'simple.in') != result['input_sha256']):
            raise ValueError('Frozen input snapshot changed during the campaign')
        result['summary'] = summarize(result['runs'])
        result['comparison_mode']='Paired committed Stage 195 versus legacy; outer-only/boundary/all-operation calibration; identical numerical work and fixed gates'
        result['status'] = 'passed'
        result.pop('active', None)
        save(result_path, result)
        print(json.dumps(result['summary'], indent=2), flush=True)
    except BaseException as error:
        result['status'] = 'failed'
        result['error'] = repr(error)
        save(result_path, result)
        raise


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'candidate', 'checked', 'legacy', 'legacy_profile', 'fixture', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=83)
    parser.add_argument('--repeats', type=int, default=6)
    parser.add_argument('--expected-checkpoints', type=int, nargs='+', default=[4])
    parser.add_argument('--require-post-restart-remap', action='store_true', default=True)
    parser.add_argument('--profiles', action='store_true', help='Retained compatibility option; paired diagnostic modes are always run')
    notes = parser.add_mutually_exclusive_group(required=True)
    notes.add_argument('--build-notes', help='Short compiler/MPI and source description')
    notes.add_argument('--build-notes-file', type=Path, help='Complete JSON build record (avoids command-line size limits)')
    return parser.parse_args(argv)


def main():
    run(parse_args())


if __name__ == '__main__':
    main()
