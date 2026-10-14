from pathlib import Path
import sys,json,hashlib,subprocess,shutil,time
R=Path('/Users/kevlahan/.codex/worktrees/connected-rhs/wavetrisk_hydrostatic');H=Path('/Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic');K=R/'.stage278-validation/rk278-connected-rhs-kit-v5'
sys.path.insert(0,str(K))
from common import runtime
from run_local_230 import validate
from analyze import analyze
from protocol_229 import ordinary
from contract import evidence
binary=Path(sys.argv[1]);out=Path(sys.argv[2]);mode=sys.argv[3];out.mkdir(parents=True)
fixture=H/'.stage276-validation/rk276-primitive-prototype-kit-v2/fixture-J4'
for n in ['simple.in','localJ4J6_checkpoint_0002.bin.zst']:shutil.copyfile(fixture/n,out/n)
(out/'grids').symlink_to(fixture/'grids',target_is_directory=True);shutil.copy2(binary,out/'climate')
env=runtime(mode=='oracle',mode in ['audit','oracle'],True);env['OMPI_MCA_mpi_yield_when_idle']='1'
if mode in ['profile','coarse']:env.update(WAVETRISK_PROFILE_BLOCK_DETAIL='1',WAVETRISK_PROFILE_PARALLEL_BLOCKS='1',WAVETRISK_PROFILE_RHS_FINE='1' if mode=='profile' else '0')
cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in'];start=time.monotonic()
with (out/'run.log').open('w') as f:subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1800)
v=validate(out,H/'.stage266-validation/local-J4/opt',4,mode=='oracle',cmd,env)
v['ordinary']=ordinary((out/'run.log').read_text(),analyze(out/'run.log'),(out/'localJ4J6_log').read_text());v['contract']=evidence(out,4,mode in ['audit','oracle'],40)
v['elapsed']=time.monotonic()-start;(out/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
print(json.dumps({'path':str(out),'status':v.get('status'),'ordinary':v['ordinary']['median'],'elapsed':v['elapsed']}))
