#!/usr/bin/env python3
"""Local J4 observer validation; timings are not cluster performance evidence."""
import argparse,json
from pathlib import Path
import local_j4
from analyze_workflow_195 import collect

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('build','fixture','out'):p.add_argument('--'+n,type=Path,required=True)
 p.add_argument('--binary',default='bin-195-opt/climate');p.add_argument('--ranks',type=int,default=4)
 p.add_argument('--mode',choices=('minimal','coarse','detail'),default='detail')
 p.add_argument('--oracles',action='store_true');a=p.parse_args();a.profile=a.mode=='detail'
 original=local_j4.environment
 def env(**kwargs):
  e=original(**kwargs);e.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',
    WAVETRISK_STAGE192_SPARSE_OUTER='0',WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',
    WAVETRISK_WORKFLOW195=str(('minimal','coarse','detail').index(a.mode)))
  return e
 local_j4.environment=env
 try:
  local_j4.run(a)
  result=collect(a.out,a.ranks,('minimal','coarse','detail').index(a.mode))
  (a.out/'workflow-summary.json').write_text(json.dumps(result,indent=2)+'\n')
 except BaseException as error:
  path=a.out/'results.json'
  if path.exists():
   r=json.loads(path.read_text());r.update(status='failed',error=repr(error));path.write_text(json.dumps(r,indent=2)+'\n')
  raise
 finally:local_j4.environment=original
if __name__=='__main__':main()
