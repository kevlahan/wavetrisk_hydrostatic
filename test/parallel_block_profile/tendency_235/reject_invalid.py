"""Inject one bounded fault into an isolated object, or poison removed scratch.

A successful scratch test must finish with every removed compatibility field
still NaN. The failure tests must hit their specific guard, not merely fail.
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
    variant='oracle' if mode=='source' else 'native'
    filename='src/parallel_block_program.f90'
    if mode=='missing':
        anchor="    if(key(3)==0) call die('missing live mass flux producer')"
        filename='src/parallel_block_diffusion.f90'
        replacement="    key=0 ! test-only dropped live producer\n"+anchor
        diagnostic='missing live mass flux producer'
    elif mode=='epoch':
        filename='src/parallel_block_velocity_owner.f90'
        # Request a different layer through the real consumer API, leaving
        # the producer's epoch guard unchanged.
        anchor='       value=owner_diffusion_value(key,k)'
        replacement='       value=owner_diffusion_value(key,k+1) ! test-only wrong epoch'
        diagnostic='external result epoch/owner'
    elif mode=='source':
        filename='src/parallel_block_program.f90'
        anchor='          o%ready=.true.'
        replacement='          if(o%code==OP_VELOCITY_SOURCE) o%value(1)=o%value(1)+1.0e-4_dp\n'+anchor
        diagnostic='independent source mismatch'
    elif mode=='scratch':
        filename='src/multi_level.f90'
        anchor='       call begin_block_diffusion(k)'
        fields=['grid(velocity_domain)%'+n+'%elts' for n in
                ('press','press_lower','geopot','geopot_lower','ke','qe','vort','bernoulli','divu')]
        fields += ['exner_fun(k)%data(velocity_domain)%elts',
                   'Laplacian_scalar(S_MASS)%data(velocity_domain)%elts',
                   'Laplacian_scalar(S_TEMP)%data(velocity_domain)%elts',
                   'Laplacian_vector(S_DIVU)%data(velocity_domain)%elts']
        poison='       do velocity_domain=1,size(grid)\n'+''.join(
            '          '+f+'=ieee_value(0.0_dp,ieee_quiet_nan)\n' for f in fields)+'       end do\n'
        replacement=poison+anchor
        diagnostic='PASS removed owner scratch remained poisoned'
    else:raise ValueError(mode)
    original=digest(tree/filename);s=(tree/filename).read_text();assert s.count(anchor)==1,(mode,anchor)
    s=s.replace(anchor,replacement)
    if mode=='scratch':
        declaration='    logical :: validate_velocity_source, overlap_mass'
        assert s.count(declaration)==1
        s=s.replace(declaration,'''    type :: Test_Primary_Backup
       real(dp),allocatable :: mass(:),temp(:),velo(:)
    end type
    type(Test_Primary_Backup),allocatable :: test_primary(:)
'''+declaration)
        start='       call begin_block_diffusion(k)'
        before='''       if(.not.allocated(test_primary)) allocate(test_primary(size(grid)))
       do velocity_domain=1,size(grid)
          test_primary(velocity_domain)%mass=q(S_MASS,k)%data(velocity_domain)%elts
          test_primary(velocity_domain)%temp=q(S_TEMP,k)%data(velocity_domain)%elts
          test_primary(velocity_domain)%velo=q(S_VELO,k)%data(velocity_domain)%elts
          q(S_MASS,k)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          q(S_TEMP,k)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
          q(S_VELO,k)%data(velocity_domain)%elts=ieee_value(0.0_dp,ieee_quiet_nan)
       end do
'''
        assert s.count(start)==1;s=s.replace(start,before+start)
        anchor='       call compute_native_velocity_gradient(q,k)'
        assert s.count(anchor)==1
        checks='\n       do velocity_domain=1,size(grid)\n'+''.join(
            "          if(.not.all(ieee_is_nan( &\n               "+f+"))) &\n               error stop 'removed owner scratch was written'\n" for f in fields)+'       end do\n'
        checks+='       do velocity_domain=1,size(grid)\n'
        for var,field in [('S_MASS','mass'),('S_TEMP','temp'),('S_VELO','velo')]:
            q='q('+var+',k)%data(velocity_domain)%elts'
            checks+="          if(.not.all(ieee_is_nan("+q+"))) error stop 'primary compatibility write'\n"
            checks+='          '+q+'=test_primary(velocity_domain)%'+field+'\n'
        checks+='       end do\n'
        s=s.replace(anchor,anchor+checks+"       if(k==zlevels) write(*,'(a)') '"+diagnostic+"'\n")
    source=out/'injected.f90';source.write_text(s)
    lines=(tree/f'build-{variant}.log').read_text().splitlines()
    compile=shlex.split(next(l for l in reversed(lines) if ' -c '+filename+' -o ' in l))
    obj=compile[compile.index('-o')+1]
    compile=[str(source) if x==filename else str(out/'injected.o') if x==obj else
             '-J'+str(out/'mod') if x.startswith('-J') else x for x in compile]
    link=shlex.split(next(l for l in reversed(lines) if obj in l and ' -c ' not in l and ' -o ' in l))
    link[link.index('-o')+1]=str(out/'climate');link=[str(out/'injected.o') if x==obj else x for x in link]
    with (out/'build.log').open('x') as f:
        for cmd in (compile,link):subprocess.run(cmd,cwd=tree,stdout=f,stderr=subprocess.STDOUT,check=True)
    text=replace_parameter((fixture/'simple.in').read_text(),'time_end','.304')
    (out/'simple.in').write_text(text);seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed)
    (out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=False,detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as f:
        code=subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
    log=(out/'run.log').read_text()
    assert (code==0 if mode=='scratch' else code!=0) and diagnostic in log,(code,log[-2000:])
    if mode=='scratch':
        from analyze import analyze
        assert analyze(out/'run.log')['completed']
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',expected_failure=mode!='scratch',returncode=code,
        diagnostic=diagnostic,original_source_sha256=original,mode=mode,command=cmd,
        primary_fields_poisoned=mode=='scratch',removed_scratch_poisoned=mode=='scratch',
        injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log')))
    print('PASS '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['missing','epoch','source','scratch'],required=True)
    run(**vars(p.parse_args()))
