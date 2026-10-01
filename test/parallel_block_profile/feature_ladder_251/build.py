"""Build a J5 dependency ladder from one sealed source and read-only main."""
import hashlib
import json
import os
import platform
from pathlib import Path
import shutil
import subprocess


EARLY = ('BLOCK_COLUMNS','BLOCK_PRIMITIVES','BLOCK_DIFFUSION',
         'BLOCK_MASS_OWNER','BLOCK_VELOCITY_OWNER','BLOCK_TENDENCY_INITIAL')
LATE = ('BLOCK_RK','PHYSICS_COLUMN_LOCAL','PHYSICS_VELOCITY_PLAN','PHYSICS_OWNER',
        'PHYSICS_OWNER_BATCH','PHYSICS_NATIVE_INPUT','PHYSICS_NATIVE_STATE',
        'SOIL_NATIVE','ADAPT_NATIVE','LIFECYCLE_NATIVE','COEFFICIENT_NATIVE',
        'MIRROR_FREE')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path,value):
    path=Path(path)
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def verify_seal(root):
    files=json.loads((root/'files-sha256.json').read_text())
    for name,expected in files.items():
        if sha(root/name)!=expected:
            raise ValueError('Sealed source differs: '+name)
    return files


def physics_objects(root):
    directory=root/'src/physics/simple_physics/phyparam/obj'
    objects={path.name:sha(path) for path in sorted(directory.glob('*.o'))}
    if len(objects)!=23:
        raise ValueError('Require 23 shared external physics objects')
    return objects


def command_for(original,stage):
    if stage not in range(7):
        raise ValueError('Unknown feature-ladder stage')
    toggles=set(EARLY)|set(LATE)
    command=[item for item in original if item.split('=',1)[0] not in toggles
             and not item.startswith(('PARAM=','BUILD_DIR=','BIN_DIR='))]
    command += [f'{name}={"on" if i<stage else "off"}' for i,name in enumerate(EARLY)]
    command += [f'{name}=off' for name in LATE]
    command += ['PARAM=param_J5',f'BUILD_DIR=build-rk251-c{stage}',f'BIN_DIR=bin-rk251-c{stage}']
    if (command[:2]!=['make','-j1'] or 'TEST_CASE=climate' not in command or
            'DEBUG=false' not in command or 'FP_CONTRACT=off' not in command or
            'STAGE202=on' not in command or 'STAGE203=on' not in command):
        raise ValueError('Candidate build selectors differ')
    return command


def run(candidate_kit,main_tree,out):
    candidate_kit,main_tree,out=map(Path,(candidate_kit,main_tree,out))
    candidate_files=verify_seal(candidate_kit)
    main_files=verify_seal(main_tree)
    identity=json.loads((main_tree/'MAIN_IDENTITY.json').read_text())
    if identity['changed_source']!=['src/checkpoint.f90']:
        raise ValueError('Unexpected main source adjustment')
    candidate=candidate_kit/'source'
    for stage in range(7):
        if (candidate/f'bin-rk251-c{stage}').exists() or (candidate/f'build-rk251-c{stage}').exists():
            raise ValueError('Use a fresh candidate source tree')
    if (main_tree/'bin-rk251-main').exists() or (main_tree/'build-rk251-main').exists():
        raise ValueError('Use a fresh main source tree')
    if any(name in os.environ for name in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES')):
        raise ValueError('Unset compiler and make overrides')
    subprocess.run(['bash',str(candidate_kit/'prepare_physics_dependencies.sh'),str(candidate)],check=True)
    original=json.loads((candidate_kit/'build-command.json').read_text())
    record=dict(host=platform.node(),main_identity=identity,
                candidate_sealed_files=len(candidate_files),main_sealed_files=len(main_files),
                compiler=subprocess.check_output(['mpif90','--version'],text=True),
                early_features=EARLY,late_features_off=LATE,builds={})

    def build_stage(stage):
        name=f'c{stage}'
        command=command_for(original,stage)
        with (out/f'build-{name}.log').open('x') as log:
            subprocess.run(command,cwd=candidate,stdout=log,stderr=subprocess.STDOUT,check=True)
        binary=candidate/f'bin-rk251-{name}/climate'
        if not binary.is_file():
            raise ValueError('Candidate stage binary absent: '+name)
        record['builds'][name]=dict(binary=str(binary),sha256=sha(binary),command=command,
                                    enabled_features=list(EARLY[:stage]))

    build_stage(0)
    shared_physics=physics_objects(candidate)
    destination=main_tree/'src/physics/simple_physics'
    if destination.is_symlink():
        destination.unlink()
    else:
        shutil.rmtree(destination)
    shutil.copytree(candidate/'src/physics/simple_physics',destination,symlinks=True)
    if physics_objects(main_tree)!=shared_physics:
        raise ValueError('Shared physics objects changed before main build')
    flags=main_tree/'rk251-flags.mk'
    flags.write_text('FFLAGS += -ffp-contract=off\n')
    main_command=['make','-f','Makefile','-f',flags.name,'-j1','TEST_CASE=climate',
                  'PARAM=param_J5','DEBUG=false','BUILD_DIR=build-rk251-main','BIN_DIR=bin-rk251-main']
    with (out/'build-main.log').open('x') as log:
        subprocess.run(main_command,cwd=main_tree,stdout=log,stderr=subprocess.STDOUT,check=True)
    main_binary=main_tree/'bin-rk251-main/climate'
    if not main_binary.is_file() or physics_objects(main_tree)!=shared_physics:
        raise ValueError('Main binary or shared physics objects differ')
    for name,expected in main_files.items():
        if name.startswith('src/physics/simple_physics/'):
            continue
        if sha(main_tree/name)!=expected:
            raise ValueError('Main non-physics source changed: '+name)
    record['builds']['main']=dict(binary=str(main_binary),sha256=sha(main_binary),
                                  command=main_command,physics_objects_sha256=shared_physics)
    save(out/'builds.json',record)
    for stage in range(1,7):
        build_stage(stage)
        if physics_objects(candidate)!=shared_physics:
            raise ValueError(f'Shared physics objects changed at feature stage c{stage}')
        save(out/'builds.json',record)
    return record
