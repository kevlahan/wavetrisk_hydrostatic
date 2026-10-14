"""MPI tests of topology-time final writers and bounded RHS phase exchange.

No fixture, solver activation, production retirement or cluster timing claim.
Runs serially, with at most four MPI ranks. A fresh output directory is required.
"""
import argparse, hashlib, json, os, subprocess
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();R=Path(__file__).resolve().parents[3];T=Path(__file__).resolve().parent
O=a.output.resolve();O.mkdir(parents=True,exist_ok=False)
files=[R/'src/kind.f90',R/'src/parallel_block_profile.f90',R/'src/parallel_block_rhs_exchange.f90',T/'exchange_unit.f90',Path(__file__).resolve(),T/'exchange_stress.f90']
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
result={'status':'building','production_active':False,'source_sha256':{str(f.relative_to(R)):sha(f) for f in files},'modes':{}}
def save(): (O/'results.json').write_text(json.dumps(result,indent=2)+'\n')
save()
failures={'missing':'missing phase producer','ambiguous':'ambiguous final writer',
 'destination':'duplicate consumer destination','orientation':'plan slot or orientation',
 'identity':'invalid producer identity or order','extent':'plan slot or orientation',
 'recompile':'release previous plan before recompiling','topology':'stale topology, mask or phase',
 'mask':'stale topology, mask or phase','phase':'stale topology, mask or phase',
 'state':'stale state or duplicate phase','duplicate':'stale state or duplicate phase',
 'distributed-mask':'inconsistent distributed plan epoch','peer-state':'peer state epoch mismatch'}
try:
 for mode in ('optimized','checked'):
  D=O/mode;D.mkdir()
  cmd=['mpifort','-ffree-line-length-132','-ffp-contract=off','-Wall','-Wextra','-Werror',
       '-J'+str(D),'-I'+str(D)]
  cmd+=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if mode=='checked' else ['-O2']
  cmd += [str(f) for f in files[:4]]+['-o',str(D/'unit')]
  print('Build exchange',mode,flush=True)
  with (D/'build.log').open('x') as f: subprocess.run(cmd,cwd=D,stdout=f,stderr=subprocess.STDOUT,check=True)
  row={'build':cmd,'binary_sha256':sha(D/'unit'),'runs':{}};result['modes'][mode]=row;save()
  for ranks,case,marker in [(1,'normal','PASS RHS exchange'),(4,'normal','PASS RHS exchange')]+[(4,k,v) for k,v in failures.items()]:
   name=f'{ranks}-{case}';cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(ranks),str(D/'unit')]
   if case!='normal':cmd.append(case)
   run=subprocess.run(cmd,cwd=D,text=True,capture_output=True,timeout=60)
   log=D/(name+'.log');log.write_text(run.stdout+run.stderr)
   ok=(run.returncode==0 if case=='normal' else run.returncode!=0) and marker in log.read_text()
   row['runs'][name]={'passed':ok,'exit_code':run.returncode,'log_sha256':sha(log)};save()
   if not ok: raise AssertionError(f'{name}: {log.read_text()}')
  stress_cmd=[str(T/'exchange_stress.f90') if x==str(T/'exchange_unit.f90') else str(D/'stress') if x==str(D/'unit') else x for x in row['build']]
  with (D/'stress-build.log').open('x') as f: subprocess.run(stress_cmd,cwd=D,stdout=f,stderr=subprocess.STDOUT,check=True)
  cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4',str(D/'stress')]
  run=subprocess.run(cmd,cwd=D,text=True,capture_output=True,timeout=60)
  log=D/'stress.log';log.write_text(run.stdout+run.stderr)
  ok=run.returncode==0 and 'PASS RHS exchange rendezvous' in log.read_text()
  row['stress']={'passed':ok,'exit_code':run.returncode,'binary_sha256':sha(D/'stress'),'log_sha256':sha(log)};save()
  if not ok: raise AssertionError(log.read_text())
 result['status']='passed';save();print('PASS optimized and checked RHS exchange tests',flush=True)
except BaseException as exc:
 result['status']='failed';result['error']=repr(exc);save();raise
