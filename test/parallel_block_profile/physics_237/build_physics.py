"""Freeze and source-bind local physics column builds. Never refresh an existing frozen tree."""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
FILES=['Makefile','src/physics/physics_simple.f90']+[str(p.relative_to(BASE)) for p in sorted((BASE/'src').glob('*.f90'))]
SELECTORS=['BLOCK_COLUMNS','BLOCK_PRIMITIVES','BLOCK_DIFFUSION','BLOCK_MASS_OWNER','BLOCK_VELOCITY_OWNER']
VARIANTS={'oracle':('check','check','on'),'native':('on','check','on'),
          'production':('on','false','on'),'control':('on','check','off')}
def build(reference,out,variants):
    out=out.resolve();reference=reference.resolve()
    if not out.exists():
        shutil.copytree(reference,out,symlinks=True,ignore=shutil.ignore_patterns('build-*','bin-*','*.log'))
        for name in FILES:
            dst=out/name
            if dst.exists() or dst.is_symlink():dst.unlink()
            shutil.copy2(BASE/name,dst,follow_symlinks=False)
    for name in FILES:
        if digest(BASE/name)!=digest(out/name):raise ValueError('Frozen source differs: '+name)
    commands={}
    for name in variants:
        owner,debug,rk=VARIANTS[name]
        cmd=['make','-j1','TEST_CASE=climate','PARAM=param_J4','DEBUG='+debug,'FP_CONTRACT=off',
             'STAGE202=on','STAGE203=on',*[s+'='+owner for s in SELECTORS],
             'BLOCK_TENDENCY_INITIAL=on','BLOCK_RK=on','PHYSICS_COLUMN_LOCAL='+('check' if name=='oracle' else rk),'BUILD_DIR=build-'+name,'BIN_DIR=bin-'+name]
        with (out/('build-'+name+'.log')).open('a') as log:
            subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
        commands[name]=cmd
    manifest=out/'build-237.json'
    previous=json.loads(manifest.read_text()) if manifest.exists() else {}
    commands={**previous.get('commands',{}),**commands}
    save(manifest,dict(stage=237,status='built',commands=commands,
         changed_sources={n:digest(out/n) for n in FILES},
         sources={str(p.relative_to(out)):digest(p) for p in sorted((out/'src').rglob('*.f90'))},
         resolved_sources={str(p.relative_to(out)):digest(p.resolve()) for p in sorted((out/'src').rglob('*.f90'))},
         binaries={n:digest(out/('bin-'+n)/'climate') for n in commands}))
    print('PASS physics frozen build '+', '.join(variants),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('reference','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=list(VARIANTS),default=list(VARIANTS))
    build(**vars(p.parse_args()))
