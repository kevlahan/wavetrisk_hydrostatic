"""Run accepted-traversal units against the current source using existing climate module interfaces."""
from pathlib import Path
import argparse,subprocess,json,hashlib
from schedule_oracle import source
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--module-build',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();R=Path(__file__).resolve().parents[3];T=Path(__file__).resolve().parent
O=a.output.resolve();O.mkdir(exist_ok=False,parents=True)
ops=subprocess.check_output(['git','show','ff740744608731a05d53771c32a3055309f3a586:src/ops.f90'],cwd=R,text=True)
(O/'oracle.f90').write_text(source(ops));result={}
cases={'':'PASS schedule','stale':'RHS slot stale layout','duplicate':'RHS schedule duplicate patch',
 'incomplete':'RHS schedule incomplete level','level':'RHS schedule level mismatch',
 'corner-order':'RHS schedule pentagon order invalid','corner-invalid':'RHS schedule corner invalid',
 'boundary-epoch':'RHS layout topology changed','missing-storage':'RHS layout missing boundary storage'}
for mode in ('optimized','checked'):
 D=O/mode;D.mkdir()
 cmd=['gfortran','-ffree-line-length-132','-ffp-contract=off','-J'+str(D),'-I'+str(D),'-I'+str(a.module_build.resolve()/'mod')]
 cmd+=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if mode=='checked' else ['-O2']
 cmd+=[str(R/'src/parallel_block_rhs_layout.f90'),str(R/'src/parallel_block_rhs_schedule.f90'),
       str(O/'oracle.f90'),str(T/'schedule_check.f90'),str(T/'schedule_unit.f90'),'-o',str(D/'unit')]
 with (D/'build.log').open('x') as f:subprocess.run(cmd,cwd=D,stdout=f,stderr=subprocess.STDOUT,check=True)
 rows={}
 for case,marker in cases.items():
  r=subprocess.run([str(D/'unit')]+([case] if case else []),cwd=D,text=True,capture_output=True)
  output=r.stdout+r.stderr;(D/((case or 'normal')+'.log')).write_text(output)
  assert (r.returncode==0 if not case else r.returncode!=0) and marker in output,output
  rows[case or 'normal']={'passed':True,'exit_code':r.returncode}
 result[mode]=rows
result['source_sha256']={str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in
 [R/'src/parallel_block_rhs_layout.f90',R/'src/parallel_block_rhs_schedule.f90',T/'schedule_check.f90',T/'schedule_unit.f90',T/'schedule_oracle.py']}
(O/'results.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS schedule units optimized and checked')
