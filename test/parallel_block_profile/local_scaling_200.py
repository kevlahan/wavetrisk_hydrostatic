#!/usr/bin/env python3
"""Local restart and observer coverage only; no cluster performance claim."""
import argparse,json
from pathlib import Path
import local_j4
from analyze_scaling_200 import collect

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','fixture','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--binary',default='bin-200-profile/climate');p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--mode',type=int,choices=(0,1),required=True)
    a=p.parse_args();a.oracles=False;a.profile=False;original=local_j4.environment
    def env(**kwargs):
        e=original(**kwargs);e.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',WAVETRISK_STAGE192_SPARSE_OUTER='0',
          WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1',
          WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200=str(a.mode));return e
    local_j4.environment=env
    try:
        local_j4.run(a);r=collect(a.out,a.ranks,a.mode)
        (a.out/'scaling-summary.json').write_text(json.dumps(r,indent=2)+'\n')
    except BaseException as error:
        q=a.out/'results.json'
        if q.exists():
            r=json.loads(q.read_text());r.update(status='failed',error=repr(error));q.write_text(json.dumps(r,indent=2)+'\n')
        raise
    finally:local_j4.environment=original
if __name__=='__main__':main()
