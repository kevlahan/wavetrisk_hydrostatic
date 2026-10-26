"""Private four-rank initialization only; never enables the256-rank candidate."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,time

def smoke(checkout,out):
    checkout=Path(checkout).resolve();out=Path(out).resolve();kit=Path(__file__).resolve().parent
    out.mkdir(parents=True);shutil.copytree(kit/'grids',out/'grids')
    text=(kit/'simple.in').read_text();lines=text.splitlines();found=0
    for i,line in enumerate(lines):
        if line.split() and line.split()[0]=='resume':lines[i]='resume -1';found+=1
    if found!=1:raise ValueError('Ambiguous fresh initialization input')
    (out/'simple.in').write_text('\n'.join(lines)+'\n')
    env={k:v for k,v in os.environ.items() if not k.startswith(('WAVETRISK_','MV2_','OMPI_','GFORTRAN_')) and not k.endswith('NUM_THREADS')}
    env.update(json.loads((kit/'checked-environment.json').read_text()))
    env.update(WAVETRISK_CAMPAIGN312_STEPS='0',WAVETRISK_SELECTIVE_BALANCE='1',WAVETRISK_SUBDIVISION_TARGET_TASKS='4',
        WAVETRISK_OWNER_MAP322='off',WAVETRISK_EXECUTION_AUDIT316='1',WAVETRISK_COMPACT_AUDIT325='1',
        WAVETRISK_WORK_LEDGER325='1',WAVETRISK_PLACEMENT322='0',WAVETRISK_WINDOW322='0',VECLIB_MAXIMUM_THREADS='1')
    (out/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('WAVETRISK_','MV2_','OMPI_','GFORTRAN_')) or k.endswith('NUM_THREADS')},indent=2)+'\n')
    binary=checkout/'bin-stage325-check/climate';t=time.monotonic()
    command=[shutil.which('mpiexec'),'-n','4',str(binary),'simple.in']
    with (out/'run.log').open('x') as log:r=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    result={'passed':r.returncode==0,'returncode':r.returncode,'seconds':time.monotonic()-t,'ranks':4,'steps':0,
            'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'native_qualification':False}
    (out/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    if r.returncode:raise RuntimeError('Initialization failed; see '+str(out/'run.log'))
    return result

if __name__=='__main__':print(json.dumps(smoke(*sys.argv[1:]),indent=2))
