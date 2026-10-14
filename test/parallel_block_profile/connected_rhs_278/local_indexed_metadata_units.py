"""Exact generic-writer reference for the preparation-only indexed metadata pull."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
files=[a.source/'src/kind.f90',a.source/'src/parallel_block_profile.f90',a.source/'src/parallel_block_rhs_exchange.f90',Path(__file__).with_name('indexed_metadata_unit.f90')]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
report=dict(status='running',sources={str(f):sha(f) for f in files},modes={})
def save():(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
try:
 for mode in ['optimized','checked']:
  d=out/mode;d.mkdir()
  cmd=['mpifort','-ffree-line-length-132','-ffp-contract=off','-Wall','-Wextra','-Werror','-J'+str(d),'-I'+str(d)]
  cmd+=['-O2'] if mode=='optimized' else ['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
  cmd+=[str(f) for f in files]+['-o',str(d/'unit')]
  with (d/'build.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
  runs=[];report['modes'][mode]=dict(build=cmd,binary_sha256=sha(d/'unit'),runs=runs)
  cases=[(1,'normal','PASS indexed metadata'),(4,'normal','PASS indexed metadata')]+[(4,k,v) for k,v in {'bad-rank':'indexed metadata home address','bad-slot':'indexed metadata source address','epoch':'inconsistent indexed metadata epoch','phase':'inconsistent indexed metadata epoch','columns':'inconsistent indexed metadata epoch'}.items()]
  for n,case,marker in cases:
   cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(n),str(d/'unit'),case]
   r=subprocess.run(cmd,capture_output=True,text=True,timeout=90);log=d/f'{n}-{case}.log';log.write_text(r.stdout+r.stderr)
   ok=(r.returncode==0 if case=='normal' else r.returncode!=0) and marker in log.read_text()
   runs.append(dict(command=cmd,case=case,passed=ok,returncode=r.returncode,log_sha256=sha(log)));save()
   if not ok:raise RuntimeError(log.read_text())
 report['status']='passed';save();print('PASS indexed metadata tests')
except BaseException as e:
 report['status']='failed';report['error']=repr(e);save();raise
