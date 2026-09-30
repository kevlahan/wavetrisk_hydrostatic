"""Reject missing/stale operands, old-field fallback and corrupt routed results.

Only an isolated replacement object/executable is built. Production files and
binaries stay unchanged. The selected guard must fire in a real MPI run.
"""
import argparse,shlex,shutil,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS

def run(tree,fixture,out,mode):
    tree,fixture,out=(p.resolve() for p in (tree,fixture,out));out.mkdir();(out/'mod').mkdir()
    variant='oracle' if mode in ('missing','corrupt') else 'production'
    if mode=='missing':
        filename='src/parallel_block_inverse.f90'
        anchor='    slot=nodes(d)%slot(id+1)\n    if(slot==0)'
        replacement='    slot=0 ! Test-only missing compiled operand\n    if(slot==0)'
        diagnostic='required scalar operand was not registered'
    elif mode=='corrupt':
        filename='src/parallel_block_columns.f90'
        anchor='          send_result(:,pos)=[lower(p),geopot(p),pressure,exner,geopot_lower]'
        replacement=anchor+'\n          send_result(4,pos)=exner*(1.0_dp+1.0e-4_dp) ! Test-only transport corruption'
        diagnostic='independent numerical mismatch: pressure layer'
    elif mode=='fallback':
        filename='src/multi_level.f90'
        anchor='       call cal_surf_press()'
        replacement='       call cal_surf_press(q) ! Test-only forbidden old-field input'
        diagnostic='block pressure production must not receive Domain solution operands'
    elif mode=='epoch':
        filename='src/parallel_block_inverse.f90'
        anchor='    provisional_operands_ready=provisional'
        replacement='    provisional_operands_ready=.false. ! Test-only stale producer epoch'
        diagnostic='provisional consumer operand epoch is not ready'
    else:raise ValueError(mode)
    original=digest(tree/filename);s=(tree/filename).read_text();assert s.count(anchor)==1
    source=out/'injected.f90';source.write_text(s.replace(anchor,replacement))
    lines=(tree/f'build-233-{variant}.log').read_text().splitlines()
    compile=shlex.split(next(l for l in reversed(lines) if ' -c '+filename+' -o ' in l))
    obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==filename else str(out/'injected.o') if x==obj else
             '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'injected.o') if x==obj else x for x in link]
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
    assert code!=0 and diagnostic in log,(code,log[-2000:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',expected_failure=True,returncode=code,
         diagnostic=diagnostic,original_source_sha256=original,mode=mode,
         injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log'),
         command=command,compile=compile,link=link))
    print('PASS '+mode+': '+diagnostic,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',required=True,choices=('missing','corrupt','fallback','epoch'))
    run(**vars(p.parse_args()))
