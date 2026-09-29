#!/usr/bin/env python3
"""Stage230 J5 production baseline, candidate and checked candidate."""
import json, os, platform, shutil, subprocess, sys, tarfile
from pathlib import Path
from capability_223 import digest, save, verify_kit, verify_sources
SPECS = {f'J5-{kind}': dict(level='J5', kind=kind) for kind in ('baseline','candidate','checked')}
KEYS = tuple(SPECS)

def source_key(kind):
    return 'reference_sources' if kind=='baseline' else 'candidate_sources'

def source_dir(kind):
    return 'reference' if kind=='baseline' else kind

def materialize(kit, root, kind):
    tree=root/source_dir(kind)
    if tree.exists():return tree
    shutil.copytree(root/'reference',tree,symlinks=True,
                    ignore=shutil.ignore_patterns('build-230-*','bin-230-*','stage230-flags.mk'))
    subprocess.run(['git','apply',str(kit/'stage230-candidate.patch')],cwd=tree,check=True)
    return tree

def build(kit, root):
    kit=kit.resolve(strict=True);root=root.resolve(strict=True);verify_kit(kit)
    m=json.loads((kit/'stage230-sources.json').read_text());path=root/'capability-230-builds.json'
    if path.exists() or any((root/source_dir(k)).exists() for k in ('candidate','checked')):
        raise ValueError('Fresh build tree required')
    if any(os.environ.get(k) for k in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES')):raise ValueError('Unset make/compiler overrides')
    verify_sources(root/'reference',m['reference_sources'])
    r=dict(stage=230,status='building',manifest_sha256=digest(kit/'stage230-sources.json'),
           build_host=platform.uname()._asdict(),builds={},fp_contract='off',
           compiler=subprocess.check_output(['mpif90','--version'],text=True),
           mpi=subprocess.run(['mpirun','--version'],capture_output=True,text=True).stdout,
           loaded_modules=os.environ.get('LOADEDMODULES'))
    save(path,r);physics=None
    try:
        subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(root/'reference')],check=True)
        for key,spec in SPECS.items():
            level,kind=spec['level'],spec['kind'];tree=materialize(kit,root,kind)
            sources=m[source_key(kind)];verify_sources(tree,sources)
            stage203='on'
            flags=tree/'stage230-flags.mk';flags.write_text('FFLAGS += -ffp-contract=off\n')
            command=['make','-f','Makefile','-f',flags.name,'-j1','TEST_CASE=climate','PARAM=param_'+level,
                     'DEBUG='+('check' if kind=='checked' else 'false'),'FP_CONTRACT=off',
                     'STAGE202=on','STAGE203='+stage203,'BUILD_DIR=build-230-'+key,'BIN_DIR=bin-230-'+key]
            log=root/f'build-230-{key}.log';r['active']=key;save(path,r)
            print('Building '+key+'; log: '+str(log),flush=True)
            with log.open('x') as f:subprocess.run(command,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
            verify_sources(tree,sources)
            current={p.name:digest(p) for p in sorted((tree/'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
            if len(current)!=23 or (physics is not None and physics!=current):raise ValueError('Physics objects differ')
            physics=current;binary=tree/f'bin-230-{key}/climate'
            r['builds'][key]=dict(binary=str(binary.relative_to(root)),sha256=digest(binary),sources=sources,
                command=command,log_sha256=digest(log),stage203=stage203,observer=kind.endswith('_observer'))
            save(path,r)
        r.update(status='passed',physics_objects=physics);r.pop('active',None);save(path,r)
        print('PASS: three fresh builds; Stage229 baseline, Stage230 candidate and checked candidate; identical physics.',flush=True)
    except BaseException as e:
        r.update(status='failed',error=repr(e));save(path,r);raise
if __name__=='__main__':build(*(Path(p) for p in sys.argv[1:]))
