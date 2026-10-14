"""Serial optimized/checked tests of private and shared exchange buffers (<=4 ranks)."""
from pathlib import Path
import argparse,subprocess,hashlib,json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();S=a.source.resolve(strict=True)/'src';T=Path(__file__).resolve().parent
D=a.output.resolve();D.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r={'status':'running','sources':{str(p):sha(p) for p in [S/'kind.f90',S/'parallel_block_profile.f90',S/'parallel_block_rhs_exchange.f90',T/'exchange_shared_unit.f90',T/'exchange_shared_stress.f90']},'runs':[]}
def save():(D/'results.json').write_text(json.dumps(r,indent=2)+'\n')
try:
 save()
 for mode in ('optimized','checked'):
  out=D/mode;out.mkdir()
  flags=['-O2'] if mode=='optimized' else ['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
  for kind in ('exchange_reference_unit','exchange_shared_unit','exchange_shared_stress','exchange_prebound_unit'):
   cmd=['mpifort','-ffree-line-length-132','-ffp-contract=off','-Wall','-Wextra','-Werror',*flags,'-J'+str(out),'-I'+str(out),str(S/'kind.f90'),str(S/'parallel_block_profile.f90'),str(S/'parallel_block_rhs_exchange.f90'),str(T/(kind+'.f90')),'-o',str(out/kind)]
   with (out/(kind+'-build.log')).open('w') as f:subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
   cases=[(4,'normal')] if kind=='exchange_shared_stress' else [(1,'normal'),(4,'normal'),(4,'peer-state'),(4,'mask'),(4,'duplicate')]
   if kind in ('exchange_shared_unit','exchange_prebound_unit'):cases.append((4,'scratch'))
   if kind=='exchange_prebound_unit':cases.extend((4,c) for c in ('provider-rank','provider-slot','mixed','binding-mode'))
   if kind in ('exchange_reference_unit','exchange_shared_unit'):cases.append((4,'ambiguous'))
   for n,case in cases:
    cmd=['mpirun','--oversubscribe','--bind-to','none','-n',str(n),str(out/kind)]+([] if case=='normal' else [case])
    run=subprocess.run(cmd,capture_output=True,text=True,timeout=60);log=out/f'{kind}-{n}-{case}.log';log.write_text(run.stdout+run.stderr)
    marker={'normal':'PASS RHS exchange','peer-state':'peer state epoch mismatch','mask':'stale topology, mask or phase','duplicate':'stale state or duplicate phase','scratch':'shared exchange scratch extent','provider-rank':'prebound provider address','provider-slot':'remote source outside workspace','mixed':'prebound provider with competing definitions','binding-mode':'inconsistent distributed plan epoch','ambiguous':'ambiguous final writer'}[case]
    ok=(run.returncode==0 if case=='normal' else run.returncode!=0) and marker in log.read_text()
    r['runs'].append(dict(mode=mode,kind=kind,ranks=n,case=case,passed=ok,log_sha256=sha(log)));save()
    if not ok:raise ValueError(log.read_text())
 r['status']='passed';save();print('PASS',len(r['runs']),'exchange cases')
except BaseException as e:r['status']='failed';r['error']=repr(e);save();raise
