#!/usr/bin/env python3
"""Build frozen Stage202, checked Stage202, matched observers and legacy.

Requires a fresh Stage196 source tree. Numerical and observer patches are
verified separately; each binary has fresh objects and identical physics.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import platform
import tarfile

VARIANTS = [('baseline','false'),('baseline_profile','false'),('checked','check')]


def digest(path):
    return hashlib.sha256(str(path.readlink()).encode() if path.is_symlink() else path.read_bytes()).hexdigest()


def verify(tree, sources):
    for name, expected in sources.items():
        if digest(tree / name) != expected:
            raise ValueError('Unexpected source content: ' + name)


def build_legacy(kit, tree, manifest, expected_physics, param="param_J5"):
    legacy = tree.parent / (tree.name+'-legacy208')
    if legacy.exists():
        raise ValueError('Fresh legacy archive directory required: '+str(legacy))
    archive=kit/'legacy-source.tar.gz'
    patch=kit/'legacy_checkpoint_io_190.patch'
    for path,key in ((archive,'legacy_archive_sha256'),(patch,'legacy_io_patch_sha256')):
        if digest(path)!=manifest[key]: raise ValueError('Legacy kit identity differs: '+path.name)
    legacy.mkdir()
    # Archive is generated from the pinned, verified repository revision.
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if not (legacy/member.name).resolve().is_relative_to(legacy.resolve()):
                raise ValueError('Unsafe legacy archive path')
            if member.issym():
                target=(legacy/member.name).parent/member.linkname
                if not target.resolve().is_relative_to(legacy.resolve()):
                    raise ValueError('Unsafe legacy archive symlink')
            if member.islnk(): raise ValueError('Unexpected hard link in legacy archive')
        tar.extractall(legacy, filter='data')
    verify(legacy,manifest['legacy_original_sources'])
    subprocess.run(['git','apply','--check',str(patch)],cwd=legacy,check=True)
    subprocess.run(['git','apply',str(patch)],cwd=legacy,check=True)
    verify(legacy,manifest['legacy_io_sources'])
    include=legacy/'src/physics/Makefile.inc'
    include.parent.mkdir(parents=True,exist_ok=True)
    if not include.exists(): shutil.copy2(kit/'physics-Makefile.inc',include)
    subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(legacy)],check=True)
    # Explicitly append the arithmetic flag: the legacy Makefile predates FP_CONTRACT.
    (legacy/'stage208-flags.mk').write_text('FFLAGS += -ffp-contract=off\n')
    command=['make','-f','Makefile','-f','stage208-flags.mk','-j1','TEST_CASE=climate','PARAM='+param,
             'DEBUG=false','BUILD_DIR=build-208-legacy','BIN_DIR=bin-208-legacy']
    log=legacy/'build-208-legacy.log'
    print('Building legacy Domain with I/O fix; log: '+str(log),flush=True)
    with log.open('x') as stream:
        subprocess.run(command,cwd=legacy,stdout=stream,stderr=subprocess.STDOUT,check=True)
    verify(legacy,manifest['legacy_io_sources'])
    physics={p.name:digest(p) for p in sorted((legacy/'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
    if physics!=expected_physics: raise ValueError('Legacy external physics differs from block builds')
    binary=legacy/'bin-208-legacy/climate'
    return {'binary':str(binary),'sha256':digest(binary),'command':command,'log_sha256':digest(log),
            'revision':manifest['legacy_revision'],'source_root':str(legacy),'param':param,
            'io_patch_sha256':digest(patch),'original_sources':manifest['legacy_original_sources'],
            'sources':manifest['legacy_io_sources'],
            'scope':'Legacy Domain numerical sources unchanged; checkpoint.f90 has only the common I/O hardening.'}


def build_legacy_profile(kit,tree,manifest,legacy_record,expected_physics):
    original=Path(legacy_record['source_root'])
    profile=tree.parent/(tree.name+'-legacy208-profile')
    if profile.exists(): raise ValueError('Fresh legacy profile directory required')
    patch=kit/'legacy-profile-208.patch'
    if digest(patch)!=manifest['legacy_profile_patch_sha256']: raise ValueError('Legacy timer patch differs')
    shutil.copytree(original,profile,symlinks=True,ignore=shutil.ignore_patterns('build-*','bin-*'))
    verify(profile,manifest['legacy_io_sources'])
    subprocess.run(['git','apply','--check',str(patch)],cwd=profile,check=True)
    subprocess.run(['git','apply',str(patch)],cwd=profile,check=True)
    verify(profile,manifest['legacy_profile_sources'])
    command=['make','-f','Makefile','-f','stage208-flags.mk','-j1','TEST_CASE=climate','PARAM=param_J5',
             'DEBUG=false','BUILD_DIR=build-208-legacy-profile','BIN_DIR=bin-208-legacy-profile']
    log=profile/'build-208-legacy-profile.log'
    print('Building instrumented legacy; log: '+str(log),flush=True)
    with log.open('x') as stream:
        subprocess.run(command,cwd=profile,stdout=stream,stderr=subprocess.STDOUT,check=True)
    verify(profile,manifest['legacy_profile_sources'])
    if {p.name:digest(p) for p in sorted((profile/'src/physics/simple_physics/phyparam/obj').glob('*.o'))}!=expected_physics:
        raise ValueError('Instrumented legacy physics differs')
    verify(original,manifest['legacy_io_sources'])
    binary=profile/'bin-208-legacy-profile/climate'
    return {'binary':str(binary),'sha256':digest(binary),'source_root':str(profile),
            'sources':manifest['legacy_profile_sources'],'command':command,'log_sha256':digest(log),
            'scope':'Experimental timer-only legacy copy; excluded from timing medians.'}


def verify_helpers(kit):
    hashes = json.loads((kit / 'helpers-sha256.json').read_text())
    for name, expected in hashes.items():
        if digest(kit / name) != expected:
            raise ValueError('Frozen helper changed: ' + name)


def build(kit, tree):
    kit, tree = kit.resolve(strict=True), tree.resolve(strict=True)
    verify_helpers(kit)
    manifest = json.loads((kit / 'stage208-sources.json').read_text())
    if any((tree.parent/(tree.name+suffix)).exists() for suffix in ('-legacy208',)):
        raise ValueError('Existing legacy208 archive; use a fresh build tree name')
    for key in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES'):
        if os.environ.get(key): raise ValueError('Unset build override before using frozen protocol: '+key)
    verify(tree, manifest['baseline_sources'])
    if digest(tree / 'Makefile') != manifest['makefile_sha256']:
        raise ValueError('Expected the validated Stage 196 Makefile')
    if digest(kit/'stage208-baseline.patch')!=manifest['reference_patch_sha256']:
        raise ValueError('Baseline patch differs')
    if digest(kit/'stage208-baseline-observer.patch')!=manifest['baseline_observer_patch_sha256']:
        raise ValueError('Baseline observer patch differs')
    if digest(kit / 'stage208-candidate.patch') != manifest['candidate_patch_sha256']:
        raise ValueError('Candidate patch hash differs')
    for kind, _ in VARIANTS:
        for prefix in ('build-208-', 'bin-208-'):
            if (tree / (prefix + kind)).exists():
                raise ValueError('Fresh build directory required: ' + str(tree / (prefix + kind)))
    path = tree / 'restriction-208-builds.json'
    if path.exists():
        raise ValueError('Existing Stage 208 build record; use a fresh checkout')
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
                'baseline_stage':202, 'reference_sources':manifest['reference_sources'],
                'baseline_sources': manifest['baseline_sources'], 'candidate_sources': manifest['candidate_sources'],
                'patch_sha256': manifest['candidate_patch_sha256'], 'makefile_sha256': manifest['makefile_sha256'],
                'compiler': compiler, 'mpi_version': mpi.stdout + mpi.stderr,
                'param': 'param_J5', 'integrator': 'RK4', 'fp_contract': 'off', 'builds': {},
                'build_host': platform.uname()._asdict(),
                'lscpu': subprocess.run(['lscpu'],capture_output=True,text=True).stdout,
                'environment': {k: os.environ[k] for k in ('FFLAGS', 'F90FLAGS', 'MPIF90', 'MAKEFLAGS', 'LOADEDMODULES')
                                if k in os.environ}}
    def save():
        path.write_text(json.dumps(identity, indent=2) + '\n')
    previous_physics = None
    save()
    try:
        for kind, debug in VARIANTS:
            if kind == 'baseline':
                patch = str(kit / 'stage208-baseline.patch')
                subprocess.run(['git', 'apply', '--check', patch], cwd=tree, check=True)
                subprocess.run(['git', 'apply', patch], cwd=tree, check=True)
            if kind=='checked':
                subprocess.run(['git','apply','--reverse',str(kit/'stage208-baseline-observer.patch')],cwd=tree,check=True)
                verify(tree,manifest['reference_sources'])
            if kind=='baseline_profile':
                subprocess.run(['git','apply',str(kit/'stage208-baseline-observer.patch')],cwd=tree,check=True)
            sources = manifest[{'baseline':'reference_sources','baseline_profile':'baseline_profile_sources','candidate':'candidate_sources','checked':'reference_sources'}[kind]]
            verify(tree, sources)
            command = ['make', '-j1', 'TEST_CASE=climate', 'PARAM=param_J5', 'DEBUG=' + debug,
                       'FP_CONTRACT=off', 'STAGE202=on', 'STAGE203=off', 'BUILD_DIR=build-208-' + kind, 'BIN_DIR=bin-208-' + kind]
            log = tree / ('build-208-' + kind.replace('_','-') + '.log')
            identity['active'] = kind
            save()
            print('Building ' + kind + '; log: ' + str(log), flush=True)
            with log.open('w') as stream:
                subprocess.run(command, cwd=tree, stdout=stream, stderr=subprocess.STDOUT, check=True)
            verify(tree, sources)
            if digest(tree / 'Makefile') != sources['Makefile']:
                raise ValueError('Makefile changed during build')
            physics = {p.name: digest(p) for p in sorted(
                (tree / 'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
            if len(physics) != 23 or (previous_physics is not None and physics != previous_physics):
                raise ValueError('Expected 23 identical external physics objects across all builds')
            previous_physics = physics
            binary = 'bin-208-' + kind + '/climate'
            identity['builds'][kind] = {'command': command, 'binary': binary, 'sha256': digest(tree / binary),
                                       'log_sha256': digest(log), 'sources':sources}
            save()
        identity['builds']['legacy'] = build_legacy(kit, tree, manifest, previous_physics)
        identity['builds']['legacy_profile']=build_legacy_profile(kit,tree,manifest,identity['builds']['legacy'],previous_physics)
        identity.update(status='passed', physics_objects=previous_physics)
        identity.pop('active', None)
        save()
        print('PASS: five fresh J5 builds: baseline, focused baseline observer, checked baseline, legacy and focused legacy observer; contraction off; identical physics objects.', flush=True)
    except BaseException as error:
        identity.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: build_restriction_208.py KIT_DIR FRESH_SOURCE_TREE')
    build(Path(sys.argv[1]), Path(sys.argv[2]))
