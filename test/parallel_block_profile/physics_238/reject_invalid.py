"""Verify the velocity plan rejects missing epochs, changed order, stale sources and corrupted operands."""
import argparse,shlex,subprocess,shutil,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS

def run(tree,fixture,out,mode):
    tree,fixture,out=(p.resolve() for p in (tree,fixture,out));out.mkdir();(out/'mod').mkdir()
    filename='src/physics/physics_simple.f90';compile_name=filename
    variant='oracle' if mode=='operand' else 'native'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if mode=='operand':
        anchor='         column_kinetic=physics_kinetic_energy(velocity_values,dom%len%elts(velocity_edges), &'
        replacement='         velocity_values(1)=velocity_values(1)+1.0_dp ! test-only corruption\n'+anchor
        diagnostic='independent physics velocity mismatch'
    elif mode=='epoch':
        anchor='    call apply_no_bdry2 (physics_call, z_null)'
        replacement='    deallocate(velocity_plan_by_domain) ! test-only missing epoch\n'+anchor
        diagnostic='physics velocity epoch absent'
    elif mode=='order':
        anchor='      if(plan%ordinal(id_i)/=consumed_columns) error stop "physics velocity consumer order mismatch"'
        replacement='      plan%ordinal(id_i)=0 ! test-only wrong order\n'+anchor
        diagnostic='physics velocity consumer order mismatch'
    elif mode=='version':
        anchor='      if(any(plan%version(plan%sources(:,id_i))/=plan%expected(:,id_i))) &'
        replacement='      plan%version(id_i)=-1 ! test-only stale source\n'+anchor
        diagnostic='physics velocity source version mismatch'
    else:raise ValueError(mode)
    assert s.count(anchor)==1;s=s.replace(anchor,replacement);source=out/'injected.f90';source.write_text(s)
    lines=(tree/f'build-{variant}.log').read_text().splitlines()
    compile=shlex.split(next(l for l in reversed(lines) if ' -c '+compile_name+' -o ' in l));obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==compile_name else str(out/'injected.o') if x==obj else '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'injected.o') if x==obj else x for x in link]
    with (out/'build.log').open('x') as f:
        for cmd in (compile,link):subprocess.run(cmd,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
    text=replace_parameter((fixture/'simple.in').read_text(),'time_end','.304');(out/'simple.in').write_text(text)
    seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed);(out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=mode=='operand',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert code!=0 and diagnostic in log,(code,log[-1600:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=True,returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=False,command=cmd))
    print('PASS physics '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['epoch','order','version','operand'],required=True);run(**vars(p.parse_args()))
