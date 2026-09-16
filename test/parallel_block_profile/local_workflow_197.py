#!/usr/bin/env python3
"""Fresh local J4/J6 restart validation; no performance inference."""
import argparse,json
from pathlib import Path
import local_j4

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('build','fixture','out'):p.add_argument('--'+key,type=Path,required=True)
 p.add_argument('--binary',default='bin-197-opt/climate');p.add_argument('--ranks',type=int,default=4)
 p.add_argument('--resident',choices=('0','1'),default='1');p.add_argument('--oracles',action='store_true')
 p.add_argument('--observe',action='store_true');a=p.parse_args();a.profile=False
 original=local_j4.environment
 def env(**kwargs):
  e=original(**kwargs);e.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',
   WAVETRISK_STAGE192_SPARSE_OUTER='0',WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',
   WAVETRISK_STAGE197_RESIDENT_WAVELETS=a.resident,WAVETRISK_WORKFLOW195='1',WAVETRISK_OPERATIONS197='2' if a.observe else '0')
  return e
 local_j4.environment=env
 try:local_j4.run(a)
 finally:local_j4.environment=original
 if a.observe:
  from analyze_operations_197 import collect
  result=collect(a.out,a.ranks,2,'block' if (a.build/'src/parallel_block_mpi.f90').exists() else 'legacy')
  (a.out/'operations-summary.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
