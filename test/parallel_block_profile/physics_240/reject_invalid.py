"""Reject invalid DAG routes/dependencies and poison old Domain consumers without restoration."""
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
        anchor='       call physics_batch_result(q(:,l),x(:,l),result(:,l))'
        replacement=anchor+'\n       result(3*n+1,l)=result(3*n+1,l)+1.0_dp ! test-only bad owner state'
        diagnostic='physics batch independent state mismatch'
    elif mode=='version':
        anchor='       do s=1,6'
        replacement='       h(9,g)=h(6,g) ! test-only future dependency\n'+anchor
        diagnostic='physics batch invalid source version'
    elif mode=='dependency':
        anchor='             do j=1,n'
        replacement='             if(round==1.and.pos==fd(sender)+1) cycle ! test-only dropped dependency\n'+anchor
        diagnostic='physics batch source epoch not ready'
    elif mode=='route':
        anchor='g=gd(rank+1)+l;r=destination(g)+1;'
        assert s.count(anchor)==2
        s=s.replace(anchor,'g=gd(rank+1)+l;r=mod(destination(g)+1,n_process)+1;')
        anchor='    deallocate(qs,xs)';replacement=anchor
        diagnostic='physics batch destination or duplicate mismatch'
    elif mode=='poison':
        s=s.replace('  use mpi_f08','  use ieee_arithmetic, only: ieee_value,ieee_quiet_nan\n  use mpi_f08',1)
        anchor='    call execute_physics_batch'
        replacement='    call apply_no_bdry2(poison_batch_sources,z_null)\n'+anchor
        extra="""  subroutine poison_batch_sources(dom,p,i,j,zlev,offs,dims,is_pole)
    type(Domain),intent(inout)::dom
    integer,intent(in)::p,i,j,zlev,offs(N_BDRY+1),dims(2,N_BDRY+1),is_pole
    integer::d,id,k
    d=dom%id+1;id=idx(i,j,offs,dims)
    do k=zmin,zlevels
       sol(S_TEMP,k)%data(d)%elts(id+1)=ieee_value(0.0_dp,ieee_quiet_nan)
       if(k>0.and.is_pole/=1) sol(S_VELO,k)%data(d)%elts(id_edge(id))=ieee_value(0.0_dp,ieee_quiet_nan)
    end do
  end subroutine poison_batch_sources
"""
        s=s.replace('  subroutine Rayleigh_friction',extra+'  subroutine Rayleigh_friction',1)
        diagnostic='Physics dependency waves rank='
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
    env=environment(oracles=mode=='output',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==(mode=='poison') and diagnostic in log,(code,log[-1600:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=mode!='poison',returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=False,command=cmd))
    print('PASS physics '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['route','dependency','version','output','poison'],required=True);run(**vars(p.parse_args()))
