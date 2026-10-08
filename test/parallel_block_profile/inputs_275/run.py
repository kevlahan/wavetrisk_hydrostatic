"""Sample actual ownership contracts with retained science gates; no timing qualification."""
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from build import build
from common import CHECKPOINTS, J5_HASHES, ORDER, VALIDATIONS, PLACEMENT_RANKS, STEP_RESTRICTIONS, preflight, runtime, save, sha, summarize
from capability_223 import coverage, numerical
from placement_275 import topology, plan, verify
from protocol_229 import ordinary
from run_local_230 import vtk_members
from validation_230 import storage_evidence
from evidence import publication_evidence
from contract import evidence as contract_evidence


def run(kit, fixture, out, repeat=False):
    if repeat: raise ValueError('No scripted repeat for diagnostic campaign')
    kit, fixture, out = map(Path, (kit, fixture, out))
    kit, fixture, out = kit.resolve(strict=True), fixture.resolve(strict=True), out.resolve()
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    identity = preflight(kit, fixture)
    if (out.exists() != repeat) or any(out.is_relative_to(p) or p.is_relative_to(out) for p in (kit, fixture)):
        raise ValueError('Use a fresh separate result directory')
    if (not os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOB_PARTITION') != 'nk' or
            int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 1 or
            socket.gethostname().split('.')[0] != 'bb13'):
        raise ValueError('Require one exclusive bb13 node through partition nk')
    for command in ('srun', 'numactl', 'lscpu', 'zstd', 'mpif90'):
        if not shutil.which(command):
            raise ValueError('Missing dependency: ' + command)
    removed = {key: os.environ.pop(key) for key in STEP_RESTRICTIONS if key in os.environ}
    out.mkdir()
    result = dict(status='running', stage=275, revision='input-contract', identity=identity, job_id=os.environ['SLURM_JOB_ID'],
                  kit_manifest_sha256=sha(kit / 'files-sha256.json'), removed_step_restrictions=removed,
                  repeat_used=False, schedule=[], runs={}, gates={}, started_utc=started_utc)
    def persist():
        result['elapsed_seconds'] = time.monotonic() - started
        save(out / 'results.json', result)
    persist()
    try:
        notes = build(kit,out)
        result['builds'] = notes
        top = topology()
        cases = plan(top)
        save(out / 'topology.json', top)
        result['host'] = top['host']
        if top['host'] != notes['host']:
            raise ValueError('Native build and execution hosts differ')
        for n in PLACEMENT_RANKS:
            save(out / f'placement-{n}.json', cases[n])
            command = launch(kit, out, n, '--probe')
            probe = subprocess.run(command, capture_output=True, text=True)
            (out / f'placement-{n}-probe.log').write_text(probe.stdout + probe.stderr)
            if probe.returncode:
                raise RuntimeError('Placement probe failed: ' + str(n))
            rows = [json.loads(line) for line in probe.stdout.splitlines() if line.startswith('{')]
            verify(rows, cases[n])
            save(out / f'placement-{n}-verified.json', rows)
        inputs = out / 'input'
        inputs.mkdir()
        for name in ('simple.in', 'test_checkpoint_0003.bin.zst'):
            shutil.copy2(fixture / name, inputs / name)
        shutil.copytree(kit / 'grids', inputs / 'grids')
        for name, expected in J5_HASHES.items():
            if sha(inputs / name) != expected:
                raise ValueError('Input snapshot changed: ' + name)
        result['input_hashes'] = J5_HASHES
        anchors = {}
        block_anchor = None

        def execute(n, kind, phase, suffix, oracle=False):
            nonlocal block_anchor
            family = 'main' if kind == 'main' else 'stage266'
            catalogue = phase=='validation' and kind in ('candidate','checked')
            label = f'J5-n{n}-{phase}-{suffix}-{kind}'
            directory = out / label
            directory.mkdir()
            for name in ('simple.in', 'test_checkpoint_0003.bin.zst'):
                shutil.copy2(inputs / name, directory / name)
            shutil.copytree(inputs / 'grids', directory / 'grids')
            record = notes['builds'][kind]
            binary = out / record['binary']
            if sha(binary) != record['binary_sha256']:
                raise ValueError('Binary changed: ' + kind)
            shutil.copy2(binary, directory / 'climate')
            result['active'] = label
            result['schedule'].append(label)
            persist()
            print('Starting ' + label, flush=True)
            command = launch(kit, out, n, directory, './climate', 'simple.in')
            with (directory / 'run.log').open('x') as log:
                code = subprocess.run(command, cwd=directory, env=runtime(oracle,catalogue=catalogue),
                                      stdout=log, stderr=subprocess.STDOUT).returncode
            row = dict(kind=kind, family=family, ranks=n, phase=phase, oracle=oracle, catalogue=catalogue,
                       returncode=code, binary_sha256=sha(directory / 'climate'),
                       input_sha256=sha(directory / 'simple.in'))
            result['runs'][label] = row
            persist()
            if code or row['binary_sha256'] != record['binary_sha256']:
                raise RuntimeError(label + ': solver/binary failure ' + str(code))
            resources = [json.loads((directory / f'rank-{rank}.json').read_text()) for rank in range(n)]
            verify(resources, cases[n])
            analysis, cov = coverage(directory, 'J5', 'test')
            if len(analysis['step_records']) != 12 or cov['remaps'] != 2:
                raise ValueError('Incomplete J5 trajectory')
            states = [step['state'] for step in analysis['step_records']]
            vtk = vtk_members(directory)
            checkpoint = sha(directory / cov['new_checkpoint'])
            if len(vtk) != 62 or checkpoint != CHECKPOINTS[family]:
                raise ValueError('Established checkpoint/VTK gate differs: ' + label)
            save(directory / 'vtk-members.json', [dict(archive=a, member=m, sha256=h)
                                                  for (a, m), h in sorted(vtk.items())])
            key = (n, family)
            if key in anchors:
                reference, old_states, old_vtk = anchors[key]
                if states != old_states or vtk != old_vtk:
                    raise ValueError('Same-rank science/VTK differs: ' + label)
                gate = numerical(reference, directory, cov['new_checkpoint'], 160, True)
            else:
                anchors[key] = (directory, states, vtk)
                reference = anchors[(16, 'main')][0]
                gate = numerical(reference, directory, cov['new_checkpoint'], 160, family == 'main')
            result['gates'][label] = gate
            save(directory / 'numerical.json', gate)
            if not gate['passed'] or states != anchors[(16, 'main')][1]:
                raise ValueError('Main-relative numerical/trajectory gate failed: ' + label)
            if family == 'stage266':
                if block_anchor is None:
                    block_anchor = directory
                block_gate = numerical(block_anchor, directory, cov['new_checkpoint'], 160, True)
                save(directory / 'block-exact.json', block_gate)
                if not block_gate['passed']:
                    raise ValueError('Cross-count exact block gate failed')
            log = (directory / 'run.log').read_text()
            row['contract']=contract_evidence(directory,n,catalogue)
            row['publication'] = publication_evidence(log, 'main' if family=='main' else 'candidate', oracle)
            if list(directory.glob('cost-rank-*.txt')):
                raise ValueError('Audit timer output in uninstrumented run')
            row['storage'] = storage_evidence(log, n, 'baseline' if family == 'main' else 'candidate', oracle)
            if family == 'stage266':
                for marker in ('Stage203 inverse workspace = on', 'Stage203 scalar fields = on'):
                    if marker not in log:
                        raise ValueError('Wrong resident inverse selectors')
                if n > 160 and f'Resident block execution: MPI tasks = {n}, geometric Domains = 160' not in log:
                    raise ValueError('Missing zero-Domain-rank evidence')
            if oracle:
                for name in ('dynamics', 'adaptation', 'vertical-remap'):
                    if f'Block {name} validation oracle enabled' not in log:
                        raise ValueError('Missing oracle: ' + name)
            if 'Parallel-block profiling summary (inclusive):' in log:
                raise ValueError('Unexpected profiling in baseline campaign')
            timing = ordinary(log, analysis, (directory / 'test_log').read_text())
            if [item['index'] for item in timing['excluded']] != [0, 7, 8]:
                raise ValueError('Ordinary-step window differs')
            row.update(ordinary=timing, checkpoint_sha256=checkpoint, steps=12, remaps=2,
                       states=states, vtk_members=62, log_sha256=sha(directory / 'run.log'),
                       vtk_manifest_sha256=sha(directory / 'vtk-members.json'))
            save(directory / 'analysis.json', dict(step_records=analysis['step_records'], coverage=cov))
            persist()
            print(f'PASS {label}: science and ownership contract; timings not qualified', flush=True)

        for index,(n,kind,oracle) in enumerate(VALIDATIONS,1):
            execute(n,kind,'validation',str(index),oracle=oracle)
        for index, (n, kind) in enumerate(ORDER,1):
            execute(n, kind, 'timing', str(index))
        rows = [result['runs'][name] for name in result['schedule'] if result['runs'][name]['phase'] == 'timing']
        result.update(status='complete', comparison=summarize(rows), ended_utc=datetime.now(timezone.utc).isoformat())
        result.pop('active', None)
        persist()
        print('COMPLETE: sampled ownership contract diagnostic; no automatic promotion.', flush=True)
    except BaseException as error:
        result.update(status='failed', error=repr(error), ended_utc=datetime.now(timezone.utc).isoformat())
        persist()
        raise


def launch(kit, out, n, *args):
    return ['srun', '--nodes=1', f'--ntasks={n}', '--cpu-bind=none', '--kill-on-bad-exit=1',
            sys.executable, str(kit / 'rank_launch_223.py'), str(out / f'placement-{n}.json'),
            *map(str, args)]


if __name__ == '__main__':
    if len(sys.argv)!=4:raise SystemExit('Usage: run.py KIT J5_FIXTURE RESULTS')
    run(*sys.argv[1:])
