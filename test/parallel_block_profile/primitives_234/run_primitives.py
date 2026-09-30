"""Source-bound owner primitive/gradient correctness test, not a performance measurement."""
import argparse,json,re,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS
from capability_223 import save,digest
from analyze import analyze

def run(binary,fixture,out,ranks=4,oracles=False,full=False,reference=None,control=False):
    binary,fixture,out=(x.resolve() for x in (binary,fixture,out))
    out.mkdir()
    text=(fixture/'simple.in').read_text()
    if not full:text=replace_parameter(text,'time_end','.304')
    (out/'simple.in').write_text(text);seed=checkpoint_name(text)
    shutil.copy2(fixture/seed,out/seed)
    (out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    shutil.copy2(binary,out/'climate')
    env=environment(oracles=oracles,detail=False);env.update(SELECTORS)
    env['OMPI_MCA_mpi_yield_when_idle']='1'
    command=['mpirun','--oversubscribe','--bind-to','none','-n',str(ranks),'./climate','simple.in']
    with (out/'run.log').open('x') as f:
        code=subprocess.run(command,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=2400).returncode
    if code:raise RuntimeError(f'{out}: solver exited {code}')
    log=(out/'run.log').read_text()
    counts=[tuple(map(int,m)) for m in re.findall(r'Block columns complete rank=(\d+) executed/checked=\s*(\d+)\s+(\d+)',log)]
    if control:
        assert not counts,counts
    else:
        assert counts and {r for r,_,_ in counts}==set(range(ranks)),counts
        assert sum(c for _,c,_ in counts)>0
    plans=[tuple(map(int,m)) for m in re.findall(r'Block columns plan rank=(\d+) geometry/owner/aliases/domains=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',log)]
    if ranks>40:assert any(owned>0 and domains==0 for r,geo,owned,aliases,domains in plans),plans
    analysis=analyze(out/'run.log');assert analysis['completed']
    if full:
        from run_local_230 import validate
        assert reference is not None
        validate(out,reference.resolve(),ranks,oracles,command,env)
    primitive=[tuple(map(int,m)) for m in re.findall(
        r'Block primitives rank=(\d+) regular/SW executed/checked=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',log)]
    pentagon=[tuple(map(int,m)) for m in re.findall(
        r'Block pentagons rank=(\d+) executed/checked=\s*(\d+)\s+(\d+)',log)]
    gradient=[tuple(map(int,m)) for m in re.findall(
        r'Block gradients rank=(\d+) executed/checked=\s*(\d+)\s+(\d+)',log)]
    gradient_plans=[tuple(map(int,m)) for m in re.findall(
        r'Block gradient plan rank=(\d+) geometry/owner=\s*(\d+)\s+(\d+)',log)]
    retention=[tuple(map(int,m)) for m in re.findall(
        r'Block gradient direct retention rank=(\d+) blocks=(\d+)',log)]
    if control:
        assert not retention
        assert not primitive and not pentagon and not gradient
    else:
        assert {r for r,n in retention}==set(range(ranks))
        for rows in (primitive,pentagon,gradient):
            assert {r[0] for r in rows}==set(range(ranks)),rows
            assert sum(r[1] for r in rows)>0
        if ranks>40:
            assert any(geo==0 and owned>0 for r,geo,owned in gradient_plans),gradient_plans
    result=dict(status='passed',control=control,
        scope='Interstage regular, southwest, pentagon and final velocity-gradient owner execution.',
        primitive_counts=primitive,pentagon_counts=pentagon,gradient_counts=gradient,gradient_plans=gradient_plans,retention=retention,
        complete_tendency_migration=False,counts=counts,plans=plans,ranks=ranks,
        command=command,binary_sha256=digest(out/'climate'),input_sha256=digest(out/'simple.in'),
        seed_sha256=digest(out/seed),log_sha256=digest(out/'run.log'),steps=len(analysis['step_records']),
        performance_claim=None)
    save(out/'primitives-results.json',result);print('PASS '+str(out),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('binary','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4);p.add_argument('--oracles',action='store_true')
    p.add_argument('--control',action='store_true')
    p.add_argument('--full',action='store_true');p.add_argument('--reference',type=Path)
    run(**vars(p.parse_args()))
