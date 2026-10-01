"""Build one changed c5 source file and compare with sealed c5 on bb13."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import socket
import statistics
import subprocess
import sys

from preflight import BASELINE_SHA, OPTIMIZED_SHA, PRIOR_JOB, SOURCE, preflight, sha

STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
PHASES = ('complete timestep', 'native velocity source', 'native velocity gradient')
STEP_RESTRICTIONS = ('SLURM_CPUS_PER_TASK', 'SRUN_CPUS_PER_TASK', 'SLURM_EXACT',
                     'SLURM_EXCLUSIVE', 'SLURM_TRES_PER_TASK', 'SLURM_HINT',
                     'SLURM_THREADS_PER_CORE')


def save(path, value):
    path = Path(path)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError('Cannot load '+str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def profile_summary(log, enabled):
    if not enabled:
        if 'Parallel-block profiling summary (inclusive):' in log:
            raise ValueError('Unprofiled run emitted a profile')
        return None
    if 'Parallel-block profiling summary (inclusive):' not in log:
        raise ValueError('Profiled run emitted no phase report')
    found = {phase: [] for phase in PHASES}
    pattern = re.compile(r'^\s*(.+?):\s*(\d+)\s*\|\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*\|')
    for line in log.splitlines():
        match = pattern.match(line)
        if match and match[1] in found:
            found[match[1]].append(dict(calls=int(match[2]), seconds_avg=float(match[4])))
    if any(len(values) != 2 for values in found.values()):
        raise ValueError('Require both complete phase-profile windows')
    return {phase: dict(windows=values, calls=sum(v['calls'] for v in values),
                        seconds_avg_sum=sum(v['seconds_avg'] for v in values))
            for phase, values in found.items()}


def run(kit, prior_kit, candidate, prior_results, out):
    kit, prior_kit, candidate, prior_results = (Path(p).resolve(strict=True)
                                               for p in (kit, prior_kit, candidate, prior_results))
    out = Path(out).resolve()
    generated_physics_changes = preflight(kit, prior_kit, candidate, prior_results)
    if out.exists() or any(out.is_relative_to(root) for root in (kit, prior_kit, candidate, prior_results)):
        raise ValueError('Use a fresh external results directory')
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 1
            or os.environ.get('SLURM_JOB_PARTITION') != 'nk' or
            socket.gethostname().split('.')[0] != 'bb13'):
        raise ValueError('Require one bb13 node through partition nk')
    for name in ('srun', 'numactl', 'lscpu', 'zstd', 'mpif90'):
        if not shutil.which(name):
            raise ValueError('Missing cluster dependency: '+name)
    if any(name in os.environ for name in ('FFLAGS', 'F90FLAGS', 'MAKEFLAGS', 'MAKEFILES')):
        raise ValueError('Unset compiler and make overrides')
    removed = {name: os.environ.pop(name) for name in STEP_RESTRICTIONS if name in os.environ}
    sys.path.insert(0, str(prior_kit))
    # Stage251 timing.py imports analyze from the source tree at import time.
    sys.path.insert(0, str(candidate/'test/parallel_block_profile'))
    prior_module = module(prior_kit/'run.py', 'rk251_runner_for_253')
    builder = module(prior_kit/'build.py', 'rk251_builder_for_253')
    from timing import ordinary
    from capability_223 import coverage
    from placement_200 import topology, verify
    from placement import plan

    prior = json.loads((prior_results/'results.json').read_text())
    top = topology()
    if top['host'] != socket.gethostname() or top['host'] != prior['host']:
        raise ValueError('Stage251 and Stage253 must use the same bb13 host')
    placement = json.loads((prior_results/'placement-83.json').read_text())
    if plan(top, 83) != placement:
        raise ValueError('Stage251 placement differs from current topology')
    verify(json.loads((prior_results/'placement-83-verified.json').read_text()), placement)
    out.mkdir()
    save(out/'topology.json', top)
    save(out/'placement-83.json', placement)
    result = dict(status='running', job_id=os.environ['SLURM_JOB_ID'], stage251_job_id=PRIOR_JOB,
                  host=top['host'], prior_results_sha256=sha(prior_results/'results.json'),
                  regenerated_physics_files=generated_physics_changes,
                  removed_step_restrictions=removed, runs={}, schedule=[], decisions={})
    save(out/'results.json', result)
    try:
        source = out/'optimized-source'
        def ignore(path, names):
            if Path(path).resolve() == candidate:
                return {name for name in names if name.startswith(('build-rk251-', 'bin-rk251-',
                                                                    'build-rk253-', 'bin-rk253-'))}
            return set()
        shutil.copytree(candidate, source, symlinks=True, ignore=ignore)
        shutil.copy2(kit/'optimized-parallel_block_program.f90', source/SOURCE)
        if sha(source/SOURCE) != OPTIMIZED_SHA:
            raise ValueError('One-file replacement failed')
        original = json.loads((candidate.parent/'build-command.json').read_text())
        command = builder.command_for(original, 5)
        command = [item.replace('build-rk251-c5', 'build-rk253-opt-c5')
                   .replace('bin-rk251-c5', 'bin-rk253-opt-c5') for item in command]
        if command[:2] != ['make', '-j1']:
            raise ValueError('Build must be serial')
        print('Building optimized c5 with make -j1', flush=True)
        with (out/'build-optimized.log').open('x') as log:
            subprocess.run(command, cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
        optimized = source/'bin-rk253-opt-c5/climate'
        if not optimized.is_file():
            raise ValueError('Optimized binary absent')
        physics_objects = builder.physics_objects(source)
        if physics_objects != prior['build']['builds']['main']['physics_objects_sha256']:
            raise ValueError('Optimized build changed Stage251 compiled physics objects')
        result['build'] = dict(command=command, source_sha256=sha(source/SOURCE),
                               baseline_binary_sha256=BASELINE_SHA,
                               optimized_binary_sha256=sha(optimized),
                               physics_objects=physics_objects)
        save(out/'results.json', result)
        env_base = prior_module.runtime_environment(candidate.parent)
        rank_launcher = candidate/'test/parallel_block_profile/rank_launch_223.py'
        input_dir = prior_results/'input'
        anchors = {}
        mapping = {item['rank']: item['cpu'] for item in placement['mapping']}

        def execute(mode, index, kind):
            label = f'J5-n83-{mode}-{index}-{kind}'
            directory = out/label
            directory.mkdir()
            for name in ('simple.in', 'test_checkpoint_0003.bin.zst'):
                shutil.copy2(input_dir/name, directory/name)
            (directory/'grids').symlink_to((input_dir/'grids').resolve(), target_is_directory=True)
            binary = (prior_results/'J5-n83-validation-c5-c5/climate') if kind == 'baseline' else optimized
            expected_binary = BASELINE_SHA if kind == 'baseline' else result['build']['optimized_binary_sha256']
            if sha(binary) != expected_binary:
                raise ValueError('Binary changed: '+kind)
            shutil.copy2(binary, directory/'climate')
            env = dict(env_base)
            env['WAVETRISK_PROFILE_PARALLEL_BLOCKS'] = '1' if mode == 'profile' else '0'
            env['WAVETRISK_PROFILE_BLOCK_DETAIL'] = '0'
            command = ['srun', '--nodes=1', '--ntasks=83', '--cpu-bind=none',
                       '--kill-on-bad-exit=1', sys.executable, str(rank_launcher),
                       str(out/'placement-83.json'), str(directory), './climate', 'simple.in']
            result['active'] = label
            result['schedule'].append(label)
            save(out/'results.json', result)
            print('Starting '+label, flush=True)
            with (directory/'run.log').open('x') as log:
                code = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                      stderr=subprocess.STDOUT).returncode
            row = dict(kind=kind, mode=mode, ranks=83, returncode=code,
                       binary_sha256=sha(directory/'climate'), log_sha256=sha(directory/'run.log'))
            result['runs'][label] = row
            save(out/'results.json', result)
            if code or row['binary_sha256'] != expected_binary:
                raise RuntimeError(f'{label} solver/binary failure {code}')
            for rank in range(83):
                record = json.loads((directory/f'rank-{rank}.json').read_text())
                if (record['rank'] != rank or record['returncode'] or record['hostname'] != top['host']
                        or record['affinity'] != [mapping[rank]]):
                    raise ValueError(f'{label} rank {rank} placement/exit mismatch')
            analysis, cov = coverage(directory, 'J5', 'test')
            if len(analysis['step_records']) != 12 or cov['remaps'] != 2:
                raise ValueError(label+': incomplete J5 trajectory')
            states = [item['state'] for item in analysis['step_records']]
            state_hash = hashlib.sha256('\n'.join(states).encode()).hexdigest()
            vtk = prior_module.vtk_members(directory)
            checkpoint = sha(directory/cov['new_checkpoint'])
            if state_hash != STATE or checkpoint != CHECKPOINT or len(vtk) != 62:
                raise ValueError(label+': J5 state/checkpoint/VTK gate differs')
            if kind in anchors and anchors[kind] != (states, vtk):
                raise ValueError(label+': repeated trajectory/VTK differs')
            anchors.setdefault(kind, (states, vtk))
            if kind == 'optimized' and vtk != anchors['baseline'][1]:
                raise ValueError(label+': optimized VTK differs from baseline')
            log = (directory/'run.log').read_text()
            prior_module.mirror_storage(log)
            row.update(steps=12, remaps=2, vtk_members=62, state_sha256=state_hash,
                       checkpoint_sha256=checkpoint,
                       ordinary=ordinary(log, analysis, (directory/'test_log').read_text()),
                       profile=profile_summary(log, mode == 'profile'))
            save(directory/'analysis.json', dict(step_records=analysis['step_records'], coverage=cov,
                                                 ordinary=row['ordinary'], profile=row['profile']))
            save(out/'results.json', result)
            print(f"PASS {label}: {row['ordinary']['median']:.9g}s", flush=True)

        for mode in ('off', 'profile'):
            for index, kind in enumerate(('baseline', 'optimized', 'optimized', 'baseline'), 1):
                execute(mode, index, kind)
            values = {kind: [row['ordinary']['median'] for row in result['runs'].values()
                             if row['mode'] == mode and row['kind'] == kind]
                      for kind in ('baseline', 'optimized')}
            result['decisions'][mode] = dict(baseline_mean_seconds=statistics.mean(values['baseline']),
                                             optimized_mean_seconds=statistics.mean(values['optimized']),
                                             difference_seconds=statistics.mean(values['optimized'])-statistics.mean(values['baseline']),
                                             runs_per_variant=2)
            save(out/'results.json', result)
        result.update(status='complete', decision='review_before_promotion')
        result.pop('active', None)
        save(out/'results.json', result)
        print('COMPLETE: one-file optimization trial; no automatic promotion.', flush=True)
    except BaseException as error:
        result.update(status='failed', error=repr(error))
        save(out/'results.json', result)
        raise


if __name__ == '__main__':
    if len(sys.argv) != 6:
        raise SystemExit('Usage: run.py KIT STAGE251_KIT CANDIDATE_SOURCE STAGE251_RESULTS FRESH_RESULTS')
    run(*sys.argv[1:])
