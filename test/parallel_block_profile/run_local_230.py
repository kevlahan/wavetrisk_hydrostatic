"""Exercise the Domain-tendency deletion through adaptation, output and reload.

Local correctness evidence only; oversubscribed timings are not performance data.
"""
import argparse
import hashlib
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

from analyze import analyze
from capability_223 import coverage, digest, numerical, save
from profile_protocol import checkpoint_name, environment, parameters
from run_local_220 import SELECTORS


def vtk_members(root):
    members = {}
    archives = sorted(root.glob('*.vtk.tgz'))
    if len(archives) != 2:
        raise ValueError('Require both output events')
    for archive in archives:
        with tarfile.open(archive) as stream:
            for member in stream:
                if member.isfile():
                    key = (archive.name, member.name)
                    if key in members:
                        raise ValueError('Duplicate output member')
                    members[key] = hashlib.sha256(stream.extractfile(member).read()).hexdigest()
    if not members:
        raise ValueError('Empty VTK outputs')
    return members


def run(binary, fixture, out, reference, ranks=4, oracles=False):
    recorded_ranks = re.findall(r'number of processors\s*=\s*(\d+)', (reference/'run.log').read_text())
    if not recorded_ranks or any(int(n) != ranks for n in recorded_ranks):
        raise ValueError('Use a completed reference at the same MPI rank count for exact VTK comparison')
    out.mkdir(exist_ok=False)
    text = (fixture/'simple.in').read_text()
    for name in ('simple.in', checkpoint_name(text)):
        shutil.copy2(fixture/name, out/name)
    (out/'grids').symlink_to((fixture/'grids').resolve(), target_is_directory=True)
    shutil.copy2(binary, out/'climate')
    env = environment(oracles=oracles, detail=False)
    env.update(SELECTORS)
    env['OMPI_MCA_mpi_yield_when_idle'] = '1'
    command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', str(ranks), './climate', 'simple.in']
    with (out/'run.log').open('x') as log:
        code = subprocess.run(command, cwd=out, env=env, stdout=log,
                              stderr=subprocess.STDOUT, timeout=1800).returncode
    if code:
        raise RuntimeError(f'{out}: exited {code}')
    return validate(out, reference, ranks, oracles, command, env)


def validate(out, reference, ranks, oracles, command, env):
    """Review a completed run; VTK requires a reference at the same rank count."""
    text = (out/'simple.in').read_text()
    p = parameters(text)
    assert (reference/'simple.in').read_bytes() == (out/'simple.in').read_bytes()
    assert digest(reference/checkpoint_name(text)) == digest(out/checkpoint_name(text))
    for directory in (out, reference):
        recorded_ranks = re.findall(r'number of processors\s*=\s*(\d+)', (directory/'run.log').read_text())
        assert recorded_ranks and all(int(n) == ranks for n in recorded_ranks), (directory, recorded_ranks)
    a, c = coverage(out, 'J4', p['run_id'])
    assert len(a['step_records']) == 12 and c['remaps'] == 2
    gate = numerical(reference, out, c['new_checkpoint'], 40, True)
    assert gate['passed'], gate
    assert [s['state'] for s in a['step_records']] == [
        s['state'] for s in analyze(reference/'run.log')['step_records']]
    log = (out/'run.log').read_text()
    absent = log.count('Domain tendency storage: absent (resident block dynamics)')
    retained = log.count('Domain tendency storage: legacy/oracle')
    assert (absent, retained) == ((0, 2) if oracles else (2, 0)), (absent, retained)
    storage_ranks = [tuple(map(int, m)) for m in
                     re.findall(r'Domain tendency storage ranks = (\d+) / (\d+)', log)]
    assert storage_ranks == [(ranks if oracles else 0, ranks)]*2, storage_ranks
    if ranks > 40:
        assert f'Resident block execution: MPI tasks = {ranks}, geometric Domains = 40' in log
    if oracles:
        for name in ('dynamics', 'adaptation', 'vertical-remap'):
            assert f'Block {name} validation oracle enabled' in log
    actual, expected = vtk_members(out), vtk_members(reference)
    assert actual == expected, 'VTK output payload differs'
    result = dict(status='passed', ranks=ranks, oracles=oracles, coverage=c,
                  reference_directory=str(reference), reference_ranks=ranks,
                  reference_binary_sha256=digest(reference/'climate'),
                  steps=len(a['step_records']), step_states_exact=True,
                  numerical=gate, vtk_exact_members=len(actual),
                  storage_headers=dict(absent=absent, retained=retained),
                  storage_ranks=storage_ranks,
                  input_sha256=digest(out/'simple.in'),
                  seed_sha256=digest(out/checkpoint_name(text)),
                  environment={k: v for k, v in env.items() if k.startswith('WAVETRISK_')},
                  command=command, binary_sha256=digest(out/'climate'), log_sha256=digest(out/'run.log'))
    save(out/'results.json', result)
    print('PASS '+str(out), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('binary', 'fixture', 'out', 'reference'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--oracles', action='store_true')
    run(**vars(parser.parse_args()))
