"""Reject stale soil epochs/ranges and corrupt aliases; poison every legacy soil operand."""
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
    filename='src/parallel_block_inverse.f90'
    readback=mode in ('owner-readback','owner-readback-41')
    success=mode=='poison' or readback
    ranks=41 if mode=='owner-readback-41' else 4
    variant='native' if mode=='alias' or readback else 'production'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if readback:
        anchor="    write(*,'(a,i0,a,i0)') 'Native soil transform rank=',rank,' layers=',inverse_levels"
        replacement="""    do c=1,2
       call check_gather203(interior_plan,c,BLOCK_PAYLOAD_SOL)
    end do
    write(*,'(a,i0)') 'Native soil owner readback exact rank=',rank
"""+anchor
        wave='       call transfer_nodes(interior_plan,c,BLOCK_PAYLOAD_WAV_COEFF,.true.)'
        assert s.count(wave)==1
        s=s.replace(wave,wave+'\n       call check_gather203(interior_plan,c,BLOCK_PAYLOAD_WAV_COEFF)')
        diagnostic='Native soil owner readback exact rank='
    elif mode=='epoch':
        anchor='    call export_native_soil_aliases(address,entry,generation)'
        replacement='    physics_results_ready=.false. ! test-only stale epoch\n'+anchor
        diagnostic='post-physics soil epoch absent'
    elif mode=='range':
        anchor='         field_first=zmin,field_count=1-zmin,entry_alias_address=address,entry_alias_value=entry)'
        replacement=anchor.replace('field_first=zmin','field_first=zmin-1')
        diagnostic='invalid inverse vertical range'
    elif mode=='alias':
        anchor="    write(*,'(a,i0,a,i0)') 'Native soil transform rank=',rank,' layers=',inverse_levels"
        replacement="""    if(size(consumer_alias_slot)>0) then
       consumer_alias_value(S_TEMP-scalars(1)+1,1,1)=consumer_alias_value(S_TEMP-scalars(1)+1,1,1)+1.0_dp
       write(*,'(a)') 'Injected soil alias corruption'
    end if
"""+anchor
        diagnostic='phase oracle mismatch'
    elif mode=='poison':
        s=s.replace('only : ieee_is_finite','only : ieee_is_finite,ieee_value,ieee_quiet_nan',1)
        s=s.replace('only : grid, sol, sol_mean','only : grid, sol, wav_coeff, sol_mean',1)
        s=s.replace('  integer :: soil_domain=0', '  type :: Soil_Backup_Type\n     real(dp),allocatable::s(:),w(:)\n  end type\n  type(Soil_Backup_Type),allocatable::soil_backup(:,:,:)\n  integer :: soil_domain=0',1)
        s=s.replace('  subroutine native_soil_transform',Path(__file__).with_name('poison.f90.inc').read_text()+'\n  subroutine native_soil_transform',1)
        anchor='    call export_native_soil_aliases(address,entry,generation)'
        replacement='    call poison_soil_fields(.false.)\n'+anchor
        marker="    write(*,'(a,i0,a,i0)') 'Native soil transform rank=',rank,' layers=',inverse_levels"
        assert s.count(marker)==1
        s=s.replace(marker,'    call poison_soil_fields(.true.)\n'+marker)
        diagnostic='Native soil Domain poison restored after transform rank='
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
    if mode!='poison':text=replace_parameter(text,'time_end','.304')
    (out/'simple.in').write_text(text)
    seed=checkpoint_name(text);shutil.copy2(fixture/seed,out/seed);(out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    env=environment(oracles=False,detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(ranks),'./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==success and diagnostic in log,(code,log[-1600:])
    if mode=='alias':assert 'Injected soil alias corruption' in log
    if readback:assert log.count(diagnostic)==ranks
    lifecycle=None
    if mode=='poison':
        assert log.count(diagnostic)==4*12
        from run_local_230 import validate
        lifecycle=validate(out,fixture,4,False,cmd,env)
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=not success,ranks=ranks,returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=(mode=='poison'),command=cmd,lifecycle=lifecycle))
    print('PASS soil '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['epoch','range','alias','poison','owner-readback','owner-readback-41'],required=True);run(**vars(p.parse_args()))
