"""Build accepted Stage265t, main and installation candidate with one shared physics build."""
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
from common import save, sha, verify


def objects(tree):
    paths = sorted((tree / 'src/physics/simple_physics/phyparam/obj').glob('*.o'))
    if len(paths) != 23:
        raise ValueError('Expected 23 external Simple-physics objects')
    return {path.name: sha(path) for path in paths}


def build(kit, out, level='J5', checked=True):
    kit, out = Path(kit).resolve(strict=True), Path(out).resolve(strict=True)
    verify(kit)
    if level not in ('J4', 'J5') or (out / 'sources').exists():
        raise ValueError('Use a fresh source directory and J4/J5')
    if any(os.environ.get(key) for key in ('FFLAGS', 'F90FLAGS', 'MAKEFLAGS', 'MAKEFILES')):
        raise ValueError('Unset make/compiler overrides')
    root = out / 'sources'
    root.mkdir()
    with tarfile.open(kit / 'sources.tar.gz') as archive:
        archive.extractall(root, filter='data')
    manifests = json.loads((kit / 'source-files.json').read_text())
    for kind in ('stage265t', 'main', 'candidate'):
        for name, expected in manifests[kind].items():
            if sha(root / kind / name) != expected:
                raise ValueError('Extracted source differs: ' + kind + '/' + name)
    notes = dict(status='building', host=platform.node(), level=level, builds={},
                 identity=json.loads((kit / 'IDENTITY.json').read_text()),
                 source_manifest_sha256=sha(kit / 'source-files.json'),
                 compiler=subprocess.check_output(['mpif90', '--version'], text=True),
                 loaded_modules=os.environ.get('LOADEDMODULES'))
    save(out / 'builds.json', notes)
    stage265t = root / 'stage265t'
    subprocess.run(['bash', str(kit / 'prepare_physics_dependencies.sh'), str(stage265t)], check=True)
    kinds = ['stage265t', 'main', 'candidate'] + (['checked'] if checked else [])
    for kind in kinds:
        tree = root / ('candidate' if kind == 'checked' else kind)
        if kind in ('main', 'candidate'):
            physics = tree / 'src/physics/simple_physics'
            shutil.rmtree(physics)
            shutil.copytree(stage265t / 'src/physics/simple_physics', physics, symlinks=True)
        flags = tree / 'stage266-flags.mk'
        flags.write_text('FFLAGS += -ffp-contract=off\n')
        command = ['make', '-f', 'Makefile', '-f', flags.name, '-j1', 'TEST_CASE=climate',
                   'PARAM=param_' + level, 'DEBUG=' + ('check' if kind == 'checked' else 'false'),
                   'BUILD_DIR=build-rk266-' + kind, 'BIN_DIR=bin-rk266-' + kind]
        if kind != 'main':
            command += ['FP_CONTRACT=off', 'STAGE202=on', 'STAGE203=on']
        notes['active'] = kind
        save(out / 'builds.json', notes)
        print('Building ' + level + '-' + kind + ' with make -j1', flush=True)
        with (out / ('build-' + kind + '.log')).open('x') as log:
            subprocess.run(command, cwd=tree, stdout=log, stderr=subprocess.STDOUT, check=True)
        binary = tree / ('bin-rk266-' + kind) / 'climate'
        current = objects(tree)
        if kind == 'stage265t':
            notes['physics_objects'] = current
        elif current != notes['physics_objects']:
            raise ValueError('Shared physics objects changed')
        source_kind = 'candidate' if kind == 'checked' else kind
        for name, expected in manifests[source_kind].items():
            if sha(tree / name) != expected:
                raise ValueError('Solver source changed during build: ' + name)
        notes['builds'][kind] = dict(command=command, binary=str(binary.relative_to(out)),
                                     binary_sha256=sha(binary), physics_objects=current,
                                     log_sha256=sha(out / ('build-' + kind + '.log')))
        save(out / 'builds.json', notes)
    notes.update(status='passed')
    notes.pop('active', None)
    save(out / 'builds.json', notes)
    return notes


def cached_build(kit, out, expected):
    """Reuse exact same solvers only after verifying sources, binaries and toolchain."""
    notes = json.loads((out / 'builds.json').read_text())
    if (notes != expected or notes['status'] != 'passed' or notes['level'] != 'J5' or
            notes['identity'] != json.loads((kit / 'IDENTITY.json').read_text()) or
            notes['source_manifest_sha256'] != sha(kit / 'source-files.json') or
            notes['loaded_modules'] != os.environ.get('LOADEDMODULES') or
            notes['compiler'] != subprocess.check_output(['mpif90', '--version'], text=True)):
        raise ValueError('Repeat build/toolchain identity differs')
    manifest = json.loads((kit / 'source-files.json').read_text())
    for kind, record in notes['builds'].items():
        source = 'candidate' if kind == 'checked' else kind
        tree = out / 'sources' / source
        if (sha(out / record['binary']) != record['binary_sha256'] or
                sha(out / ('build-' + kind + '.log')) != record['log_sha256'] or
                objects(tree) != notes['physics_objects']):
            raise ValueError('Repeat binary/log/physics changed: ' + kind)
        for name, digest in manifest[source].items():
            if sha(tree / name) != digest:
                raise ValueError('Repeat source changed: ' + name)
    print('PASS: optional repeat reuses four verified binaries; no rebuild', flush=True)
    return notes
