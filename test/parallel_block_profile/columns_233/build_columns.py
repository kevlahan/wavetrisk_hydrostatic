"""Build the bounded column-owner integration in a fresh frozen tree."""
from pathlib import Path
import argparse,shutil,subprocess,sys
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
FILES=['Makefile','src/parallel_block_columns.f90','src/parallel_block_primitive_kernels.f90',
       'src/parallel_block_inverse.f90','src/parallel_block_mpi.f90','src/parallel_block.f90',
       'src/diagnostics.f90','src/multi_level.f90','src/time_integr.f90']
def build(reference,out):
    reference,out=reference.resolve(),out.resolve()
    if out.exists():
        for name in FILES:
            if digest(BASE/name)!=digest(out/name):raise ValueError('Frozen source differs: '+name)
    else:
        shutil.copytree(reference,out,symlinks=True)
        for name in FILES:shutil.copy2(BASE/name,out/name)
    commands=[]
    for name,selector,debug in [('oracle','check','check'),('production','on','false'),('control','off','check')]:
        cmd=['make','-j1','TEST_CASE=climate','PARAM=param_J4','DEBUG='+debug,'FP_CONTRACT=off',
             'STAGE202=on','STAGE203=on','BLOCK_COLUMNS='+selector,
             'BUILD_DIR=build-233-'+name,'BIN_DIR=bin-233-'+name]
        with (out/('build-233-'+name+'.log')).open('a') as f:
            code=subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT).returncode
        if code:raise RuntimeError('build failed: '+name)
        commands.append(cmd)
    sources={str(p.relative_to(out)):digest(p) for p in (out/'src').rglob('*.f90')}
    save(out/'build-233.json',dict(status='passed',scope='Interstage column owner integration',
        commands=commands,changed_sources={n:digest(out/n) for n in FILES},sources=sources,
        binaries={n:digest(out/('bin-233-'+n)/'climate') for n in ('oracle','production','control')}))
    print('PASS three frozen builds',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('reference','out'):p.add_argument('--'+n,type=Path,required=True)
    build(**vars(p.parse_args()))
