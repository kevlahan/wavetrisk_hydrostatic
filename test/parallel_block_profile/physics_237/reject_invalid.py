"""Verify removal of physics hydrostatic scratch without saving/restoring an image."""
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
    variant='oracle' if mode=='source' else 'native'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if mode=='source':
        anchor='      column_surface=grav_accel*column_surface+p_top'
        replacement=anchor+'\n      column_surface=column_surface+1.0_dp ! test-only corruption'
        diagnostic='independent physics hydrostatic mismatch'
    elif mode=='epoch':
        anchor='    call apply_no_bdry2 (physics_call, z_null)'
        replacement='    deallocate(column_visits) ! test-only missing epoch\n'+anchor
        diagnostic='Simple-physics column epoch absent'
    elif mode=='duplicate':
        anchor='      if(column_visits(d)%seen(id_i)/=0) error stop "duplicate Simple-physics column consumer"'
        replacement='      column_visits(d)%seen(id_i)=1 ! test-only repeated consumer\n'+anchor
        diagnostic='duplicate Simple-physics column consumer'
    elif mode=='poison':
        s=s.replace('  implicit none','  use ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan\n  implicit none',1)
        s=s.replace('    integer :: d,p','    integer :: d,p,k',1)
        fields=['grid(d)%'+x+'%elts' for x in ('surf_press','press_lower','press','geopot','geopot_lower')]
        anchor='    call apply_no_bdry2 (physics_call, z_null)'
        poison='    do d=1,size(grid)\n'+''.join('       '+f+'=ieee_value(0.0_dp,ieee_quiet_nan)\n' for f in fields)
        poison+='       do k=1,zlevels\n          exner_fun(k)%data(d)%elts=ieee_value(0.0_dp,ieee_quiet_nan)\n       end do\n    end do\n'
        check='\n    do d=1,size(grid)\n'+''.join("       if(.not.all(ieee_is_nan("+f+"))) error stop 'removed physics scratch written'\n" for f in fields)
        check+="       do k=1,zlevels\n          if(.not.all(ieee_is_nan(exner_fun(k)%data(d)%elts))) error stop 'removed physics Exner written'\n       end do\n    end do\n"
        diagnostic='PASS physics scratch remained poisoned'
        replacement=poison+anchor+check+"    write(*,'(a)') '"+diagnostic+"'\n"
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
    env=environment(oracles=mode=='source',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==(mode=='poison') and diagnostic in log,(code,log[-1600:])
    if mode=='poison':
        from analyze import analyze
        assert analyze(out/'run.log')['completed'] and log.count(diagnostic)==4
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=mode!='poison',returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=False,command=cmd))
    print('PASS physics '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['poison','source','epoch','duplicate'],required=True);run(**vars(p.parse_args()))
