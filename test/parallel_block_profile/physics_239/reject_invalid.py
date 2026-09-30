"""Reject invalid owner requests and prove package soil is reloaded per column."""
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
    variant='oracle' if mode=='output' else 'native'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if mode=='output':
        anchor='    call execute_owner_packet(header,reference)'
        replacement='    packet(4*zlevels+2)=packet(4*zlevels+2)+1.0_sp ! test-only output corruption\n'+anchor
        diagnostic='independent physics owner output mismatch'
    elif mode=='version':
        anchor='    if(any(header(9:14)<0).or.any(header(9:14)>=header(6))) error stop "physics owner invalid source version"'
        # Corrupt at producer, not inside the receiver validation routine.
        anchor='      header(9:14)=plan%expected(:,id)'
        replacement=anchor+'\n      header(9)=header(6) ! test-only future source version'
        diagnostic='physics owner invalid source version'
    elif mode=='order':
        anchor='    if(header(6)<=last_owner_column(source+1)) error stop "physics owner duplicate or reordered request"'
        replacement='    last_owner_column(source+1)=header(6) ! test-only consumed request\n'+anchor
        diagnostic='physics owner duplicate or reordered request'
    elif mode=='route':
        anchor='    if(b>0) destination=block_catalog(b)%owner'
        replacement=anchor+'\n    if(b>0) destination=mod(destination+1,n_process) ! test-only wrong destination'
        diagnostic='physics owner destination mismatch'
    elif mode=='soil':
        s=s.replace('  use mpi_f08','  use ieee_arithmetic, only: ieee_value,ieee_quiet_nan\n  use soil_mod, only: package_surface=>Tsurf,package_soil=>Tsoil\n  use mpi_f08',1)
        anchor='    first=header(7)==1'
        replacement=anchor+"\n    if(allocated(package_surface)) package_surface=ieee_value(0.0_sp,ieee_quiet_nan)\n    if(allocated(package_soil)) package_soil=ieee_value(0.0_sp,ieee_quiet_nan)\n    if(.not.first) owner_package_started=.false. ! test-only late first-use path\n"
        diagnostic='Physics owner rank='
    else:raise ValueError(mode)
    assert s.count(anchor)==1;s=s.replace(anchor,replacement);source=out/'injected.f90';source.write_text(s)
    lines=(tree/f'build-{variant}.log').read_text().splitlines()
    compile=shlex.split(next(l for l in reversed(lines) if ' -c '+compile_name+' -o ' in l));obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==compile_name else str(out/'injected.o') if x==obj else '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'injected.o') if x==obj else x for x in link]
    with (out/'build.log').open('x') as f:
        for cmd in (compile,link):subprocess.run(cmd,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
    text=replace_parameter((fixture/'simple.in').read_text(),'time_end','.3065' if mode=='soil' else '.304');(out/'simple.in').write_text(text)
    seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed);(out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=mode=='output',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==(mode=='soil') and diagnostic in log,(code,log[-1600:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=mode!='soil',returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=False,command=cmd))
    print('PASS physics '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['route','order','version','output','soil'],required=True);run(**vars(p.parse_args()))
