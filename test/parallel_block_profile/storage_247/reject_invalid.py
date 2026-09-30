"""Reject stale lifecycle/coefficient epochs and accidental persistent legacy views."""
import argparse,shlex,subprocess,shutil,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS
from analyze import analyze

def run(tree,fixture,out,mode):
    tree,fixture,out=(p.resolve() for p in (tree,fixture,out));out.mkdir();(out/'mod').mkdir()
    success=False;ranks=4;variant='production'
    filename='src/domain.f90' if mode=='mirror' else 'src/parallel_block_inverse.f90'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if mode=='lifecycle':
        anchor="    if(.not.native_lifecycle_ready().or.generation/=plan_generation) call die('initial native alias generation')"
        replacement='    lifecycle_ready=.false. ! injected stale epoch\n'+anchor
        diagnostic='initial native alias generation'
    elif mode=='coefficients':
        anchor="    if(.not.native_coefficients_ready) call die('native coefficient handoff absent')"
        replacement='    native_coefficients_ready=.false. ! injected missing handoff\n'+anchor
        diagnostic='native coefficient handoff absent'
    elif mode=='mirror':
        anchor="    if(domain_prognostic_view_present()) error stop 'persistent Domain prognostic storage present'"
        replacement='    call materialize_domain_prognostic_view ! injected forbidden mirror\n'+anchor
        diagnostic='persistent Domain prognostic storage present'
    else:raise ValueError(mode)
    assert s.count(anchor)==1;s=s.replace(anchor,replacement)
    source=out/'injected.f90';source.write_text(s)
    lines=(tree/f'build-{variant}.log').read_text().splitlines()
    compile=shlex.split(next(l for l in reversed(lines) if ' -c '+filename+' -o ' in l));obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==filename else str(out/'injected.o') if x==obj else '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'injected.o') if x==obj else x for x in link]
    with (out/'build.log').open('x') as f:
        for cmd in (compile,link):subprocess.run(cmd,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
    text=(fixture/'simple.in').read_text()
    text=replace_parameter(text,'time_end','.304')
    (out/'simple.in').write_text(text)
    seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed);(out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=False,detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(ranks),'./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==success and diagnostic in log,(code,log[-1600:])
    lifecycle=None
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=not success,ranks=ranks,returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=False,command=cmd,lifecycle=lifecycle))
    print('PASS storage guard '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['lifecycle','coefficients','mirror'],required=True);run(**vars(p.parse_args()))
