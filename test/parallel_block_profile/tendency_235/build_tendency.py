"""Freeze and build the complete climate tendency kernel migration.

The supplied reference must already contain the climate physics library.
Existing frozen trees are verified, never silently refreshed.
"""
import argparse,json,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
SELECTORS=['BLOCK_COLUMNS','BLOCK_PRIMITIVES','BLOCK_DIFFUSION','BLOCK_MASS_OWNER','BLOCK_VELOCITY_OWNER']
FILES=['Makefile']+[str(p.relative_to(BASE)) for p in sorted((BASE/'src').glob('*.f90'))]
VARIANTS={'oracle':('check','check'),'native':('on','check'),'production':('on','false'),'control':('off','check')}
def build(reference,out,variants):
    reference,out=reference.resolve(),out.resolve()
    if not out.exists():
        shutil.copytree(reference,out,symlinks=True,ignore=shutil.ignore_patterns('build-*','bin-*','*.log'))
        for name in FILES:shutil.copy2(BASE/name,out/name)
    for name in FILES:
        if digest(BASE/name)!=digest(out/name):raise ValueError('Frozen source differs: '+name)
    commands={}
    for name in variants:
        selector,debug=VARIANTS[name]
        cmd=['make','-j1','TEST_CASE=climate','PARAM=param_J4','DEBUG='+debug,'FP_CONTRACT=off',
             'STAGE202=on','STAGE203=on',*[s+'='+selector for s in SELECTORS],
             'BLOCK_TENDENCY_INITIAL='+('off' if selector=='off' else 'on'),
             'BUILD_DIR=build-'+name,'BIN_DIR=bin-'+name]
        with (out/('build-'+name+'.log')).open('a') as f:
            subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
        commands[name]=cmd
    manifest=out/'build-235.json'
    prior=json.loads(manifest.read_text()) if manifest.exists() else {}
    commands={**prior.get('commands',{}),**commands}
    save(manifest,dict(status='passed',commands=commands,
        sources={str(p.relative_to(out)):digest(p) for p in sorted((out/'src').rglob('*.f90'))},
        source_contents={str(p.relative_to(out)):digest(p.resolve()) for p in sorted((out/'src').rglob('*.f90'))},
        changed_sources={n:digest(out/n) for n in FILES},
        binaries={n:digest(out/('bin-'+n)/'climate') for n in commands}))
    print('PASS frozen builds: '+', '.join(variants),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('reference','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--variants',nargs='+',choices=list(VARIANTS),default=list(VARIANTS))
    build(**vars(p.parse_args()))
