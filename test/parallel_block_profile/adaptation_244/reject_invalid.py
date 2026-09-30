"""Reject stale adaptation epochs/generations and verify native owners without Domain operands."""
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
    readback=mode in ('owner-readback','owner-readback-41')
    success=mode=='poison' or readback
    ranks=41 if mode=='owner-readback-41' else 4
    variant='native' if readback else 'production'
    filename='src/main.f90' if mode=='poison' else 'src/parallel_block_inverse.f90'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if readback:
        anchor="    write(*,'(a,i0)') 'Native adaptation reconstruction rank=',rank"
        replacement="""    do c=1,2
       call check_gather203(interior_plan,c,BLOCK_PAYLOAD_SOL)
       call check_gather203(interior_plan,c,BLOCK_PAYLOAD_WAV_COEFF)
    end do
    write(*,'(a,i0)') 'Native adaptation owner readback exact rank=',rank
"""+anchor
        diagnostic='Native adaptation owner readback exact rank='
    elif mode=='epoch':
        anchor="    if(.not.native_transition_ready()) call die('native adaptation mask entry')"
        replacement='    transition_phase=0 ! injected stale epoch\n'+anchor
        diagnostic='native adaptation mask entry'
    elif mode=='generation':
        anchor="    if(generation<=transition_generation) call die('native adaptation generation did not advance')"
        replacement='    transition_generation=generation ! injected stale generation\n'+anchor
        diagnostic='native adaptation generation did not advance'
    elif mode=='poison':
        s=s.replace('type(Initial_State), allocatable :: ini_st(:)',"""type(Initial_State), allocatable :: ini_st(:)
  type :: Adapt_Backup_Type
     real(dp),allocatable::s(:),w(:)
  end type
  type(Adapt_Backup_Type),allocatable::adapt_backup(:,:,:)
""",1)
        s=s.replace('contains','contains\n'+Path(__file__).with_name('poison.f90.inc').read_text(),1)
        anchor='    call prepare_native_transition_mask'
        replacement='    call poison_adaptation_fields(.false.)\n'+anchor
        marker='       call native_adaptation_reconstruct(native_velocity_plan_generation())'
        assert s.count(marker)==1
        s=s.replace(marker,marker+'\n       call poison_adaptation_fields(.true.)')
        diagnostic='Native adaptation poison restored rank='
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
    print('PASS adaptation '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['epoch','generation','poison','owner-readback','owner-readback-41'],required=True);run(**vars(p.parse_args()))
