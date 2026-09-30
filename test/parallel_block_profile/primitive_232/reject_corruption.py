"""Prove the primitive comparison detects an arithmetic corruption.

Compile one test-only kernel object and link against the checked comparison
build. Neither the original kernel source nor its executable is overwritten.
"""
import argparse,shlex,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS

def run(tree,fixture,out):
    tree,fixture,out=(p.resolve() for p in (tree,fixture,out));out.mkdir();(out/'mod').mkdir()
    filename='src/parallel_block_primitive_kernels.f90';original=digest(tree/filename)
    s=(tree/filename).read_text();anchor='    exner=cp*(pressure/p0)**kappa'
    assert s.count(anchor)==1
    s=s.replace(anchor,anchor+'\n    exner=exner*(1.0_dp+1.0e-4_dp) ! Deliberate test-only corruption')
    source=out/'corrupted.f90';source.write_text(s)
    lines=(tree/'build-232.log').read_text().splitlines()
    compile=shlex.split(next(l for l in lines if ' -c '+filename+' -o ' in l))
    obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==filename else str(out/'corrupted.o') if x==obj else
             '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'corrupted.o') if x==obj else x for x in link]
    with (out/'build.log').open('x') as f:
        for command in (compile,link):subprocess.run(command,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
    text=replace_parameter((fixture/'simple.in').read_text(),'time_end','.304')
    (out/'simple.in').write_text(text);seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed)
    (out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=False,detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    command=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as f:
        code=subprocess.run(command,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600).returncode
    log=(out/'run.log').read_text()
    assert code!=0 and 'primitive232 numerical mismatch' in log and 'column pressure' in log,(code,log[-2000:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',expected_failure=True,returncode=code,
         diagnostic='primitive232 numerical mismatch: column pressure',original_source_sha256=original,
         injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log'),
         command=command,compile=compile,link=link))
    print('PASS: corrupted pressure kernel rejected',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    run(**vars(p.parse_args()))
