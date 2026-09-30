"""One real adaptive-fixture step for all primitive kernel comparisons."""
import argparse,re,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS
from capability_223 import save,digest
from analyze import analyze

def run(tree,fixture,out):
    tree,fixture,out=(x.resolve() for x in (tree,fixture,out))
    out.mkdir()
    text=replace_parameter((fixture/'simple.in').read_text(),'time_end','.304')
    (out/'simple.in').write_text(text);seed=checkpoint_name(text)
    shutil.copy2(fixture/seed,out/seed)
    (out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    shutil.copy2(tree/'bin-230-checked/climate',out/'climate')
    env=environment(oracles=False,detail=False);env.update(SELECTORS)
    env['OMPI_MCA_mpi_yield_when_idle']='1'
    command=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as f:
        code=subprocess.run(command,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600).returncode
    log=(out/'run.log').read_text()
    if code:raise RuntimeError(f'{out}: solver exited {code}')
    counts=[[int(v) for v in line.split()] for line in re.findall(r'Primitive232 exact counts:([^\n]*)',log)]
    assert len(counts)==4 and all(len(c)==10 for c in counts),counts
    totals=[sum(c[i] for c in counts) for i in range(10)]
    assert all(v>0 for v in totals),totals
    assert len(analyze(out/'run.log')['step_records'])==1
    result=dict(status='passed',scope='One-step arithmetic extraction check; geometry-owner execution retained.',
        counts=dict(zip(['regular','southwest','pentagon_NE','pentagon_SE','pentagon_SW','pentagon_NW',
                        'pressure_layer','surface_pressure','scalar_diffusion','velocity_diffusion'],totals)),
        command=command,binary_sha256=digest(out/'climate'),input_sha256=digest(out/'simple.in'),
        seed_sha256=digest(out/seed),log_sha256=digest(out/'run.log'),
        production_owner_migration=False,performance_claim=None)
    save(out/'results.json',result);print(result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    run(**vars(p.parse_args()))
