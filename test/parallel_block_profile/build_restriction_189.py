#!/usr/bin/env python3
"""Build Stage 188 baseline then Stage 189 candidate in a fresh isolated checkout.

The supplied tree is intentionally patched only after the baseline build passes.
No clean, reset, commit, fetch, submission or changes to another checkout occur.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

VARIANTS = [('baseline', 'false'), ('candidate', 'false'), ('checked', 'check')]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(tree, sources):
    for name, expected in sources.items():
        if digest(tree / name) != expected:
            raise ValueError('Unexpected source content: ' + name)


def build(kit, tree):
    kit, tree = kit.resolve(strict=True), tree.resolve(strict=True)
    manifest = json.loads((kit / 'stage189-sources.json').read_text())
    verify(tree, manifest['baseline_sources'])
    if digest(tree / 'Makefile') != manifest['makefile_sha256']:
        raise ValueError('Expected the validated Stage 188 Makefile')
    if digest(kit / 'stage189-candidate.patch') != manifest['candidate_patch_sha256']:
        raise ValueError('Candidate patch hash differs')
    for kind, _ in VARIANTS:
        for prefix in ('build-189-', 'bin-189-'):
            if (tree / (prefix + kind)).exists():
                raise ValueError('Fresh build directory required: ' + str(tree / (prefix + kind)))
    path = tree / 'restriction-189-builds.json'
    if path.exists():
        raise ValueError('Existing Stage 189 build record; use a fresh checkout')
    include = tree / 'src/physics/Makefile.inc'
    if include.exists() and digest(include) != digest(kit / 'physics-Makefile.inc'):
        raise ValueError('Unexpected physics interface Makefile.inc')
    include.parent.mkdir(parents=True, exist_ok=True)
    if not include.exists():
        shutil.copy2(kit / 'physics-Makefile.inc', include)
    subprocess.run(['bash', str(kit / 'prepare_physics_dependencies.sh'), str(tree)], check=True)
    compiler = subprocess.check_output(['mpif90', '--version'], text=True)
    mpi = subprocess.run(['mpirun', '--version'], capture_output=True, text=True)
    identity = {'status': 'building', 'baseline_commit': manifest['baseline_commit'],
                'baseline_sources': manifest['baseline_sources'], 'candidate_sources': manifest['candidate_sources'],
                'patch_sha256': manifest['candidate_patch_sha256'], 'makefile_sha256': manifest['makefile_sha256'],
                'compiler': compiler, 'mpi_version': mpi.stdout + mpi.stderr,
                'param': 'param_J5', 'integrator': 'RK4', 'fp_contract': 'off', 'builds': {},
                'environment': {k: os.environ[k] for k in ('FFLAGS', 'F90FLAGS', 'MPIF90', 'MAKEFLAGS', 'LOADEDMODULES')
                                if k in os.environ}}
    def save():
        path.write_text(json.dumps(identity, indent=2) + '\n')
    previous_physics = None
    save()
    try:
        for kind, debug in VARIANTS:
            if kind == 'candidate':
                patch = str(kit / 'stage189-candidate.patch')
                subprocess.run(['git', 'apply', '--check', patch], cwd=tree, check=True)
                subprocess.run(['git', 'apply', patch], cwd=tree, check=True)
            sources = manifest['baseline_sources'] if kind == 'baseline' else manifest['candidate_sources']
            verify(tree, sources)
            command = ['make', '-j1', 'TEST_CASE=climate', 'PARAM=param_J5', 'DEBUG=' + debug,
                       'FP_CONTRACT=off', 'BUILD_DIR=build-189-' + kind, 'BIN_DIR=bin-189-' + kind]
            log = tree / ('build-189-' + kind + '.log')
            identity['active'] = kind
            save()
            print('Building ' + kind + '; log: ' + str(log), flush=True)
            with log.open('w') as stream:
                subprocess.run(command, cwd=tree, stdout=stream, stderr=subprocess.STDOUT, check=True)
            verify(tree, sources)
            if digest(tree / 'Makefile') != manifest['makefile_sha256']:
                raise ValueError('Makefile changed during build')
            physics = {p.name: digest(p) for p in sorted(
                (tree / 'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
            if len(physics) != 23 or (previous_physics is not None and physics != previous_physics):
                raise ValueError('Expected 23 identical external physics objects across all builds')
            previous_physics = physics
            binary = 'bin-189-' + kind + '/climate'
            identity['builds'][kind] = {'command': command, 'binary': binary, 'sha256': digest(tree / binary),
                                       'log_sha256': digest(log)}
            save()
        identity.update(status='passed', physics_objects=previous_physics)
        identity.pop('active', None)
        save()
        print('PASS: three fresh J5 builds, contraction off; identical external physics objects.', flush=True)
    except BaseException as error:
        identity.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: build_restriction_189.py KIT_DIR FRESH_SOURCE_TREE')
    build(Path(sys.argv[1]), Path(sys.argv[2]))
