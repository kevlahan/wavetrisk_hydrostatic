#!/usr/bin/env python3
"""Build three fresh Stage 187 arithmetic controls from a verified source tree."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(kit, tree):
    kit, tree = kit.resolve(strict=True), tree.resolve(strict=True)
    manifest = json.loads((kit / 'arithmetic-188-sources.json').read_text())
    for name, expected in manifest['sources'].items():
        if digest(tree / name) != expected:
            raise ValueError('Source differs from validated Stage 187: ' + name)
    variants = [('baseline', 'default', 'false', 'default'),
                ('candidate', 'off', 'false', 'off'),
                ('checked', 'check-off', 'check', 'off')]
    for _, label, _, _ in variants:
        for prefix in ('build-188-', 'bin-188-'):
            if (tree / (prefix + label)).exists():
                raise ValueError('Fresh build directory required; already exists: ' + str(tree / (prefix + label)))
    makefile = tree / 'Makefile'
    supported = {manifest['makefile_after']}
    if manifest.get('makefile_promoted'):
        supported.add(manifest['makefile_promoted'])
    actual = digest(makefile)
    if actual == manifest['makefile_before']:
        patch = kit / 'fp_contract_188.patch'
        if digest(patch) != manifest['patch_sha256']:
            raise ValueError('Build-option patch hash differs')
        subprocess.run(['git', 'apply', '--check', str(patch)], cwd=tree, check=True)
        subprocess.run(['git', 'apply', str(patch)], cwd=tree, check=True)
    elif actual not in supported:
        raise ValueError('Makefile differs from the supported validated versions')
    if digest(makefile) not in supported:
        raise ValueError('Patched Makefile hash differs')
    include = tree / 'src/physics/Makefile.inc'
    if not include.exists():
        shutil.copy2(kit / 'physics-Makefile.inc', include)
    subprocess.run(['bash', str(kit / 'prepare_physics_dependencies.sh'), str(tree)], check=True)
    identity = {'sources': manifest['sources'], 'makefile_sha256': digest(makefile),
                'compiler': subprocess.check_output(['mpif90', '--version'], text=True),
                'param': 'param_J5', 'integrator': 'RK4', 'builds': {}}
    previous_physics = None
    for kind, label, debug, contract in variants:
        command = ['make', '-j1', 'PARAM=param_J5', 'DEBUG=' + debug,
                   'FP_CONTRACT=' + contract, 'BUILD_DIR=build-188-' + label,
                   'BIN_DIR=bin-188-' + label]
        log = tree / ('build-188-' + label + '.log')
        print('Building ' + label + '; log: ' + str(log), flush=True)
        with log.open('w') as stream:
            subprocess.run(command, cwd=tree, stdout=stream, stderr=subprocess.STDOUT, check=True)
        physics = {p.name: digest(p) for p in sorted(
            (tree / 'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
        if not physics or (previous_physics is not None and physics != previous_physics):
            raise ValueError('External physics objects changed between comparison builds')
        previous_physics = physics
        binary = 'bin-188-' + label + '/climate'
        identity['builds'][kind] = {'command': command, 'binary': binary,
                                   'sha256': digest(tree / binary), 'log_sha256': digest(log)}
    for name, expected in manifest['sources'].items():
        if digest(tree / name) != expected:
            raise ValueError('Source changed during build: ' + name)
    identity['physics_objects'] = previous_physics
    (tree / 'arithmetic-188-builds.json').write_text(json.dumps(identity, indent=2) + '\n')
    print('PASS: three fresh builds; identical external physics objects.', flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: build_arithmetic_188.py KIT_DIR SOURCE_TREE')
    build(Path(sys.argv[1]), Path(sys.argv[2]))
