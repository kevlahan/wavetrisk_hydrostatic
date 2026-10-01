"""Build sealed migrated and read-only main sources on the same bb13 host."""
import hashlib
import json
import platform
from pathlib import Path
import shutil
import subprocess


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def verify_seal(root):
    files = json.loads((root/'files-sha256.json').read_text())
    for name, expected in files.items():
        if sha(root/name) != expected:
            raise ValueError('Sealed source differs: '+name)
    return files


def physics_objects(root):
    directory = root/'src/physics/simple_physics/phyparam/obj'
    objects = {path.name:sha(path) for path in sorted(directory.glob('*.o'))}
    if len(objects) != 23:
        raise ValueError('Require 23 shared external physics objects')
    return objects


def run(candidate_kit, main_tree, out, perf):
    candidate_kit, main_tree, out, perf = map(Path, (candidate_kit, main_tree, out, perf))
    candidate_files = verify_seal(candidate_kit)
    main_files = verify_seal(main_tree)
    identity = json.loads((main_tree/'MAIN_IDENTITY.json').read_text())
    if identity['changed_source'] != ['src/checkpoint.f90']:
        raise ValueError('Unexpected main source adjustment')
    candidate = candidate_kit/'source'
    for tree in (candidate, main_tree):
        if (tree/'bin-rk250').exists() or (tree/'build-rk250').exists():
            raise ValueError('Use fresh unbuilt source trees')
    if any(k in __import__('os').environ for k in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES')):
        raise ValueError('Unset compiler and make overrides')
    patch = perf/'quiet_storage.patch'
    patch_identity = json.loads((perf/'PATCH.json').read_text())
    if sha(patch) != patch_identity['patch_sha256']:
        raise ValueError('Storage diagnostic patch differs')
    if candidate_files['source/src/domain.f90'] != patch_identity['original_domain_sha256']:
        raise ValueError('Frozen candidate Domain source differs')
    subprocess.run(['git','apply','--unidiff-zero','--check',str(patch)],cwd=candidate,check=True)
    subprocess.run(['git','apply','--unidiff-zero',str(patch)],cwd=candidate,check=True)
    if sha(candidate/'src/domain.f90') != patch_identity['patched_domain_sha256']:
        raise ValueError('Patched candidate Domain source differs')
    subprocess.run(['bash', str(candidate_kit/'prepare_physics_dependencies.sh'), str(candidate)], check=True)
    commands = json.loads((candidate_kit/'build-command.json').read_text())
    commands = [x for x in commands if not x.startswith(('PARAM=','BUILD_DIR=','BIN_DIR='))]
    commands += ['PARAM=param_J5','BUILD_DIR=build-rk250','BIN_DIR=bin-rk250']
    if commands[:2] != ['make','-j1'] or 'MIRROR_FREE=on' not in commands or 'FP_CONTRACT=off' not in commands:
        raise ValueError('Candidate build selectors differ')
    record = dict(host=platform.node(), main_identity=identity,
                  candidate_sealed_files=len(candidate_files), main_sealed_files=len(main_files),
                  diagnostic_patch=patch_identity,
                  compiler=subprocess.check_output(['mpif90','--version'], text=True), builds={})
    with (out/'build-candidate.log').open('x') as log:
        subprocess.run(commands, cwd=candidate, stdout=log, stderr=subprocess.STDOUT, check=True)
    candidate_binary = candidate/'bin-rk250/climate'
    if not candidate_binary.is_file():
        raise ValueError('Candidate binary absent')
    candidate_physics = physics_objects(candidate)
    record['builds']['candidate'] = dict(binary=str(candidate_binary),sha256=sha(candidate_binary),command=commands)

    # Stage229's direct-main protocol used one external Simple-physics build for
    # both solvers. Keep the numerical solver sources separate and share only it.
    destination = main_tree/'src/physics/simple_physics'
    if destination.is_symlink():
        destination.unlink()
    else:
        shutil.rmtree(destination)
    shutil.copytree(candidate/'src/physics/simple_physics', destination, symlinks=True)
    if physics_objects(main_tree) != candidate_physics:
        raise ValueError('Shared physics objects changed before main build')
    flags = main_tree/'rk250-flags.mk'
    flags.write_text('FFLAGS += -ffp-contract=off\n')
    main_command = ['make','-f','Makefile','-f',flags.name,'-j1','TEST_CASE=climate',
                    'PARAM=param_J5','DEBUG=false','BUILD_DIR=build-rk250','BIN_DIR=bin-rk250']
    with (out/'build-main.log').open('x') as log:
        subprocess.run(main_command, cwd=main_tree, stdout=log, stderr=subprocess.STDOUT, check=True)
    main_binary = main_tree/'bin-rk250/climate'
    if not main_binary.is_file() or physics_objects(main_tree) != candidate_physics:
        raise ValueError('Main binary or shared physics identity differs')
    for name, expected in main_files.items():
        if name.startswith('src/physics/simple_physics/'):
            continue
        if sha(main_tree/name) != expected:
            raise ValueError('Main non-physics source changed: '+name)
    record['builds']['main'] = dict(binary=str(main_binary),sha256=sha(main_binary),command=main_command,
                                    physics_objects_sha256=candidate_physics)
    save(out/'builds.json', record)
    return record
