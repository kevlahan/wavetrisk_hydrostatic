"""Isolated RK producer faults and destructive old-primary poison, without backups."""
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
    filename='src/parallel_block_inverse.f90'
    if mode=='missing':
        anchor="    if(.not.rk_base_ready.or..not.rk_mass_ready(k)) call die('RK mass alias layer unavailable')"
        replacement='    rk_mass_ready(k)=.false. ! test-only missing producer\n'+anchor
        diagnostic='RK mass alias layer unavailable'
    elif mode=='epoch':
        anchor="    if(rk_boundary_ready.or.stage/=rk_stage+1) call die('RK closure stage sequence')"
        replacement='    rk_stage=rk_stage+1 ! test-only stale epoch\n'+anchor
        diagnostic='RK closure stage sequence'
    elif mode=='source':
        anchor='          transfer_recv(q+1:q+EDGE)=vector_closure(interior_plan%key(:,p),k,h)'
        replacement=anchor+'\n          transfer_recv(q+1)=transfer_recv(q+1)+1.0_dp ! test-only corrupt producer'
        diagnostic='independent RK alias mismatch'
    elif mode=='base':
        filename='src/parallel_block.f90'
        anchor='  if(.not.block_tendency_commit_checkpoint_ready) error stop "RK base checkpoint unavailable"'
        replacement='  block_tendency_commit_checkpoint_ready=.false. ! test-only missing base\n'+anchor
        diagnostic='RK base checkpoint unavailable'
    elif mode=='fallback':
        filename='src/time_integr.f90'
        anchor='          if(validate_oracle) call RK_sub_step_compatibility(q,h/4,q1)'
        replacement='          call RK_sub_step_compatibility(q,h/4,q1) ! test-only forbidden fallback'
        diagnostic='removed Domain RK compatibility closure requested'
    elif mode=='poison':
        anchor='    rk_stage=stage;rk_stages=stages'
        replacement='''    ! Test-only: destroy original registered atmospheric inputs after binding.
    ! There is deliberately no saved Domain image and no restoration.
    do p=1,original_inverse_nodes
       d=node_domain(p);id=node_id(p)
       do k=1,zlevels
          do n=scalars(1),scalars(2)
             if(stage>1) then
                if(.not.ieee_is_nan(sol(n,k)%data(d)%elts(id+1))) call die('old RK scalar was written')
             end if
             sol(n,k)%data(d)%elts(id+1)=ieee_value(0.0_dp,ieee_quiet_nan)
          end do
          if(stage>1) then
             if(.not.all(ieee_is_nan(sol(S_VELO,k)%data(d)%elts(EDGE*id+1:EDGE*id+EDGE)))) &
                  call die('old RK vector was written')
          end if
          sol(S_VELO,k)%data(d)%elts(EDGE*id+1:EDGE*id+EDGE)=ieee_value(0.0_dp,ieee_quiet_nan)
       end do
    end do
    if(stage==stages) write(*,'(a,i0)') 'PASS old RK primary stayed poisoned rank=',rank
'''+anchor
        diagnostic='PASS old RK primary stayed poisoned rank='
    else:raise ValueError(mode)
    original=digest(tree/filename);s=(tree/filename).read_text();assert s.count(anchor)==1,(mode,anchor)
    s=s.replace(anchor,replacement)
    if mode=='poison':s=s.replace('use ieee_arithmetic, only : ieee_is_finite',
                                  'use ieee_arithmetic, only : ieee_is_finite,ieee_is_nan,ieee_value,ieee_quiet_nan')
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
    env=environment(oracles=mode=='source',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as f:
        code=subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
    log=(out/'run.log').read_text()
    assert (code==0 if mode=='poison' else code!=0) and diagnostic in log,(code,log[-2000:])
    if mode=='poison':
        from analyze import analyze
        assert analyze(out/'run.log')['completed']
        assert log.count(diagnostic)==4
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',expected_failure=mode!='poison',returncode=code,
        diagnostic=diagnostic,original_source_sha256=original,mode=mode,command=cmd,
        primary_fields_poisoned=mode=='poison',domain_backup=False,
        injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),log_sha256=digest(out/'run.log')))
    print('PASS '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['missing','epoch','source','base','fallback','poison'],required=True)
    run(**vars(p.parse_args()))
