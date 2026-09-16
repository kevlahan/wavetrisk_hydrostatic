#!/usr/bin/env python3
"""Local correctness/coverage run of unchanged numerical code with Stage198 probes."""
import argparse,json
from pathlib import Path
import local_j4

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('build','fixture','out'):p.add_argument('--'+key,type=Path,required=True)
 p.add_argument('--binary',required=True);p.add_argument('--ranks',type=int,default=4)
 p.add_argument('--mode',type=int,choices=range(4),default=1);p.add_argument('--offset',type=int,default=0)
 a=p.parse_args();a.profile=False;a.oracles=False
 original=local_j4.environment
 def env(**kwargs):
  e=original(**kwargs);e.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',WAVETRISK_STAGE192_SPARSE_OUTER='0',WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1',WAVETRISK_WORKFLOW195='1',WAVETRISK_OPERATIONS198=str(a.mode),WAVETRISK_SAMPLE_OFFSET198=str(a.offset));return e
 local_j4.environment=env
 try:
  local_j4.run(a)
  from analyze_operations_198 import collect
  from analyze_workflow_195 import collect as workflow
  side='block' if (a.build/'src/parallel_block_mpi.f90').exists() else 'legacy'
  r=collect(a.out,a.ranks,a.mode,side);workflow(a.out,a.ranks,1)
  (a.out/'operations-summary.json').write_text(json.dumps(r,indent=2)+'\n')
 except BaseException as error:
  q=a.out/'results.json'
  if q.exists():
   r=json.loads(q.read_text());r.update(status='failed',error=repr(error));q.write_text(json.dumps(r,indent=2)+'\n')
  raise
 finally:local_j4.environment=original
if __name__=='__main__':main()
