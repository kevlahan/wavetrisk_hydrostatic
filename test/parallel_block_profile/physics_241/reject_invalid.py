"""Reject invalid native input routes/epochs and poison old input arrays during gathering."""
import argparse,shlex,subprocess,shutil,sys,re
from pathlib import Path
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'test/parallel_block_profile'))
from capability_223 import digest,save
from profile_protocol import environment,replace_parameter,checkpoint_name
from run_local_220 import SELECTORS

def run(tree,fixture,out,mode):
    tree,fixture,out=(p.resolve() for p in (tree,fixture,out));out.mkdir();(out/'mod').mkdir()
    filename='src/physics/physics_simple.f90';compile_name=filename
    variant='oracle' if mode=='native-value' else 'native'
    s=(tree/filename).read_text();original=digest(tree/filename)
    if mode=='native-value':
        anchor='    call MPI_Alltoallv(owned,nv*rc,nv*rd,MPI_DOUBLE_PRECISION,received,nv*sc,nv*sd, &'
        replacement='    owned(1,:)=owned(1,:)+10.0_dp ! test-only corrupt scalar inputs\n    owned(4*n+1:7*n,:)=owned(4*n+1:7*n,:)+0.125_dp ! corrupt edge inputs\n'+anchor
        diagnostic='physics batch independent packet mismatch'
    elif mode=='native-route':
        anchor='    call compile_local_block_inverse_routes(recvkey,address)'
        replacement="""    if(size(recvkey,2)>0) then
       do b=1,size(block_catalog)
          if(block_catalog(b)%owner/=rank) exit
       end do
       recvkey(1,1)=b ! test-only wrong native storage owner
    end if
"""+anchor
        diagnostic='inverse routes: nonlocal block'
    elif mode=='native-epoch':
        anchor='    batch_phase=1'
        replacement='    call finish_physics_consumer_epoch ! test-only stale epoch\n'+anchor
        diagnostic='physics consumer alias epoch is not ready'
    elif mode=='entry-poison':
        s=s.replace('  use mpi_f08','  use ieee_arithmetic, only: ieee_value,ieee_quiet_nan\n  use mpi_f08',1)
        s=s.replace('  integer,allocatable :: batch_entry_key(:,:,:)','  type :: Entry_Backup_Type\n     real(dp),allocatable :: mass(:,:),temp(:,:),velo(:,:),mean_m(:,:),mean_t(:,:),soil(:,:)\n  end type\n  type(Entry_Backup_Type),allocatable :: entry_backup(:)\n  integer,allocatable :: batch_entry_key(:,:,:)',1)
        s=s.replace('  subroutine fill_native_physics_entry',Path(__file__).with_name('entry_poison.f90.inc').read_text()+'\n  subroutine fill_native_physics_entry',1)
        s=s.replace('    call execute_physics_batch','    call execute_physics_batch\n    call restore_native_entry',1)
        anchor='    batch_phase=1'
        replacement=anchor+'\n    call poison_native_entry'
        diagnostic='Native physics entry poison restored after gather'
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
    env=environment(oracles=mode=='native-value',detail=False);env.update(SELECTORS);env['OMPI_MCA_mpi_yield_when_idle']='1'
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
    with (out/'run.log').open('x') as log:code=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
    log=(out/'run.log').read_text();assert (code==0)==(mode=='entry-poison') and diagnostic in log,(code,log[-1600:])
    assert digest(tree/filename)==original
    save(out/'results.json',dict(status='passed',mode=mode,expected_failure=mode!='entry-poison',returncode=code,diagnostic=diagnostic,
         original_source_sha256=original,injected_source_sha256=digest(source),binary_sha256=digest(out/'climate'),
         log_sha256=digest(out/'run.log'),input_sha256=digest(out/'simple.in'),seed_sha256=digest(out/seed),domain_backup=(mode=='entry-poison'),command=cmd,
         primary_alias_columns=[list(map(int,m)) for m in re.findall(r'Native physics primary aliases rank=(\d+) columns=(\d+)',log)]))
    print('PASS physics '+mode,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('tree','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--mode',choices=['native-route','native-epoch','native-value','entry-poison'],required=True);run(**vars(p.parse_args()))
