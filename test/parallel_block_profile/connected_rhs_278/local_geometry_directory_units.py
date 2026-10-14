"""Use matching complete solver object builds to test real geometry-directory types."""
from pathlib import Path
import argparse,subprocess,shlex,json,hashlib
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--optimized-log',type=Path,required=True)
p.add_argument('--checked-log',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
unit=Path(__file__).with_name('geometry_directory_unit.f90');sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
report={'status':'running','unit_sha256':sha(unit),'modes':{}}
def save():(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
try:
 for mode,log,build in [('optimized',a.optimized_log,'build-probe'),('checked',a.checked_log,'build-structural-check')]:
  d=out/mode;d.mkdir()
  lines=[s for s in log.read_text().splitlines() if s.startswith('mpif90 ') and ' -o bin-' in s]
  cmd=shlex.split(lines[-1]);i=cmd.index('-o');cmd[i+1]=str(d/'unit');cmd=[x for x in cmd if x!=build+'/test.o']
  cmd[1:1]=['-ffree-line-length-132','-ffp-contract=off','-I'+str(a.source/build/'mod'),'-J'+str(d),str(unit)]
  cmd[1:1]=['-O2'] if mode=='optimized' else ['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
  with (d/'build.log').open('w') as f:subprocess.run(cmd,cwd=a.source,stdout=f,stderr=subprocess.STDOUT,check=True)
  row={'build':cmd,'binary_sha256':sha(d/'unit'),'runs':[]};report['modes'][mode]=row
  cases=[(1,'normal','PASS geometry directory'),(4,'normal','PASS geometry directory')]+[(4,k,v) for k,v in {'source-key':'RHS geometry directory source identity','request-key':'RHS geometry request identity','layout':'RHS changed geometry directory layout','stale':'RHS stale geometry directory'}.items()]
  for n,case,marker in cases:
   command=['mpirun','--oversubscribe','--bind-to','none','-n',str(n),str(d/'unit'),case]
   run=subprocess.run(command,capture_output=True,text=True,timeout=90);f=d/f'{n}-{case}.log';f.write_text(run.stdout+run.stderr)
   ok=(run.returncode==0 if case=='normal' else run.returncode!=0) and marker in f.read_text()
   row['runs'].append(dict(command=command,passed=ok,returncode=run.returncode,log_sha256=sha(f)));save()
   if not ok:raise RuntimeError(f.read_text())
 report['status']='passed';save();print('PASS geometry directory tests')
except BaseException as e:
 report['status']='failed';report['error']=repr(e);save();raise
