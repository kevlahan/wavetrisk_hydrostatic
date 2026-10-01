"""Build matched mirror and mirror-free binaries serially on the allocated nk host."""
import hashlib
import json
import platform
from pathlib import Path
import subprocess


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def verify_seal(kit):
    files = json.loads((kit/'files-sha256.json').read_text())
    for name, expected in files.items():
        if sha(kit/name) != expected:
            raise ValueError('Candidate kit file differs: '+name)
    return len(files)


def commands(kit, level, comparison):
    if level not in ('J4', 'J5') or comparison not in ('mirrored','feature-off'):
        raise ValueError((level,comparison))
    original = json.loads((kit/'build-command.json').read_text())
    result = {}
    for kind, mirror in (('baseline', 'off'), ('candidate', 'on')):
        command = [x for x in original if not x.startswith(('PARAM=', 'BUILD_DIR=', 'BIN_DIR=', 'MIRROR_FREE='))]
        command += ['PARAM=param_'+level, 'MIRROR_FREE='+mirror,
                    f'BUILD_DIR=build-{kind}-{level}', f'BIN_DIR=bin-{kind}-{level}']
        if kind=='baseline' and comparison=='feature-off':
            disabled={'BLOCK_RK','PHYSICS_COLUMN_LOCAL','PHYSICS_VELOCITY_PLAN','PHYSICS_OWNER',
                      'PHYSICS_OWNER_BATCH','PHYSICS_NATIVE_INPUT','PHYSICS_NATIVE_STATE',
                      'SOIL_NATIVE','ADAPT_NATIVE','LIFECYCLE_NATIVE','COEFFICIENT_NATIVE'}
            command=[(arg.partition('=')[0]+'=off') if arg.partition('=')[0] in disabled else arg for arg in command]
        assert command[:2] == ['make', '-j1'] and 'FP_CONTRACT=off' in command and 'DEBUG=false' in command
        result[kind] = command
    return result


def run(kit, out, level, comparison):
    kit, out = Path(kit).resolve(), Path(out).resolve()
    count = verify_seal(kit)
    source = kit/'source'
    commands_by_kind = commands(kit, level, comparison)
    subprocess.run(['bash', str(kit/'prepare_physics_dependencies.sh'), str(source)], check=True)
    versions = {}
    for name, args in [('gfortran', ['gfortran', '--version']),
                       ('mpi', ['mpif90', '--version']), ('netcdf', ['nf-config', '--all'])]:
        versions[name] = subprocess.check_output(args, text=True)
    builds = {}
    for kind, command in commands_by_kind.items():
        binary = source/f'bin-{kind}-{level}/climate'
        if binary.exists():
            raise ValueError('A build already exists: '+str(binary))
        with (out/f'build-{kind}.log').open('x') as log:
            subprocess.run(command, cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
        if not binary.is_file():
            raise ValueError('Missing built binary: '+str(binary))
        builds[kind] = dict(command=command, binary=str(binary), sha256=sha(binary))
    result = dict(level=level, comparison=comparison, host=platform.node(), sealed_file_count=count,
                  local_validation_sha256=sha(kit/'LOCAL_VALIDATION.json'),
                  candidate_kit_manifest_sha256=sha(kit/'files-sha256.json'),
                  versions=versions, builds=builds)
    save(out/'builds.json', result)
    return result
