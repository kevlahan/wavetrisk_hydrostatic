#!/usr/bin/env python3
"""Build three fresh binaries: uninstrumented A/B and route-audited checked B."""
import json,os,platform,subprocess,sys,shutil
from pathlib import Path
from build_restriction_208 import digest,verify,verify_helpers

KINDS=('baseline','candidate','checked')

def build(kit,tree):
    kit=kit.resolve(strict=True);tree=tree.resolve(strict=True);verify_helpers(kit)
    m=json.loads((kit/'stage215-sources.json').read_text());verify(tree,m['baseline_sources'])
    if any(os.environ.get(k) for k in ('FFLAGS','F90FLAGS','MAKEFLAGS','MAKEFILES')):
        raise ValueError('Unset compiler/make overrides')
    path=tree/'restriction-215-builds.json'
    if path.exists() or any((tree/(p+k)).exists() for k in KINDS for p in ('build-215-','bin-215-')):
        raise ValueError('Fresh tree/build directories required')
    include=tree/'src/physics/Makefile.inc';include.parent.mkdir(parents=True,exist_ok=True)
    if include.exists() and digest(include)!=digest(kit/'physics-Makefile.inc'):raise ValueError('Physics include differs')
    if not include.exists():shutil.copy2(kit/'physics-Makefile.inc',include)
    subprocess.run(['bash',str(kit/'prepare_physics_dependencies.sh'),str(tree)],check=True)
    (tree/'stage215-flags.mk').write_text('FFLAGS += -ffp-contract=off\n')
    (tree/'stage215-audit.mk').write_text('FFLAGS += -DSTAGE215_AUDIT\n')
    r=dict(stage=215,status='building',fp_contract='off',builds={},param='param_J5',integrator='RK4',
           build_host=platform.uname()._asdict(),compiler=subprocess.check_output(['mpif90','--version'],text=True),
           manifest_sha256=digest(kit/'stage215-sources.json'),
           mpi=subprocess.run(['mpirun','--version'],capture_output=True,text=True).stdout,
           loaded_modules=os.environ.get('LOADEDMODULES'))
    def save():path.write_text(json.dumps(r,indent=2)+'\n')
    save();physics=None
    try:
        for kind in KINDS:
            if kind in ('baseline','candidate'):
                patch=kit/('stage215-baseline.patch' if kind=='baseline' else 'stage215-candidate.patch')
                subprocess.run(['git','apply','--check',str(patch)],cwd=tree,check=True)
                subprocess.run(['git','apply',str(patch)],cwd=tree,check=True)
            sources=m['reference_sources' if kind=='baseline' else 'candidate_sources'];verify(tree,sources)
            command=['make','-f','Makefile','-f','stage215-flags.mk']
            if kind=='checked':command+=['-f','stage215-audit.mk']
            command+=['-j1','TEST_CASE=climate','PARAM=param_J5','DEBUG='+('check' if kind=='checked' else 'false'),
                      'FP_CONTRACT=off','STAGE202=on','STAGE203=off','BUILD_DIR=build-215-'+kind,'BIN_DIR=bin-215-'+kind]
            log=tree/('build-215-'+kind+'.log');r['active']=kind;save();print('Building '+kind+'; '+str(log),flush=True)
            with log.open('x') as f:subprocess.run(command,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
            verify(tree,sources)
            current={p.name:digest(p) for p in sorted((tree/'src/physics/simple_physics/phyparam/obj').glob('*.o'))}
            if len(current)!=23 or (physics is not None and current!=physics):raise ValueError('Physics objects differ')
            physics=current;binary='bin-215-'+kind+'/climate'
            r['builds'][kind]=dict(binary=binary,sha256=digest(tree/binary),command=command,sources=sources,log_sha256=digest(log));save()
        r.update(status='passed',physics_objects=physics);r.pop('active',None);save()
        print('PASS: three fresh J5 builds; uninstrumented baseline/candidate and route-audited checked candidate; identical physics.',flush=True)
    except BaseException as e:r.update(status='failed',error=repr(e));save();raise

if __name__=='__main__':build(*(Path(p) for p in sys.argv[1:]))
