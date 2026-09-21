#!/usr/bin/env python3
"""Six fresh builds on one host: J4/J5 Stage222 baseline, candidate and checked candidate."""
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from capability_223 import digest, save, verify_kit, verify_sources
from candidate_capability_223 import candidate

KEYS = tuple(f'{level}-{kind}' for level in ('J4','J5') for kind in ('baseline','candidate','checked'))


def build(kit, root):
    kit=kit.resolve(strict=True);root=root.resolve(strict=True);verify_kit(kit)
    m=json.loads((kit/'stage223-sources.json').read_text())
    path=root/'capability-223-builds.json';reference=root/'reference';trial=root/'candidate'
    if path.exists() or trial.exists():raise ValueError('Fresh build tree required')
    if any(os.environ.get(k) for k in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES')):
        raise ValueError('Unset compiler/make overrides')
    verify_sources(reference,m['reference_sources'])
    r=dict(stage=223,status='building',manifest_sha256=digest(kit/'stage223-sources.json'),
           build_host=platform.uname()._asdict(),builds={},fp_contract='off',
           compiler=subprocess.check_output(['mpif90','--version'],text=True),
           mpi=subprocess.run(['mpirun','--version'],capture_output=True,text=True).stdout,
           loaded_modules=os.environ.get('LOADEDMODULES'))
    save(path,r);physics=None
    try:
        subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(reference)],check=True)
        for key in KEYS:
            level,kind=key.split('-');tree=reference if kind=='baseline' else trial
            if kind!='baseline' and not trial.exists():
                # Carry the already-built external physics objects into the
                # candidate; every solver still gets fresh object directories.
                shutil.copytree(reference,trial,symlinks=True,
                                ignore=shutil.ignore_patterns('build-223-*','bin-223-*'))
                candidate(trial);verify_sources(trial,m['candidate_sources'])
            sources=m['reference_sources' if kind=='baseline' else 'candidate_sources']
            verify_sources(tree,sources)
            flags=tree/'stage223-flags.mk';flags.write_text('FFLAGS += -ffp-contract=off\n')
            command=['make','-f','Makefile','-f',flags.name,'-j1','TEST_CASE=climate','PARAM=param_'+level,
                     'DEBUG='+('check' if kind=='checked' else 'false'),'FP_CONTRACT=off',
                     'STAGE202=on','STAGE203=off','BUILD_DIR=build-223-'+key,'BIN_DIR=bin-223-'+key]
            log=root/('build-223-'+key+'.log');r['active']=key;save(path,r)
            print('Building '+key+'; log: '+str(log),flush=True)
            with log.open('x') as f:subprocess.run(command,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
            verify_sources(tree,sources)
            current={p.name:digest(p) for p in sorted((tree/'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
            if len(current)!=23 or (physics is not None and physics!=current):raise ValueError('Physics objects differ')
            physics=current;binary=tree/('bin-223-'+key)/'climate'
            r['builds'][key]=dict(binary=str(binary.relative_to(root)),sha256=digest(binary),sources=sources,
                                  command=command,log_sha256=digest(log));save(path,r)
        r.update(status='passed',physics_objects=physics);r.pop('active',None);save(path,r)
        print('PASS: six fresh J4/J5 builds, checked candidates, contraction off, identical physics.',flush=True)
    except BaseException as e:
        r.update(status='failed',error=repr(e));save(path,r);raise


if __name__=='__main__':build(*(Path(p) for p in sys.argv[1:]))
