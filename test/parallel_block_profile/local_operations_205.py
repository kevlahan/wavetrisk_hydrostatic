#!/usr/bin/env python3
"""Four-rank J4 restart validation with the frozen Stage201 algorithm switches."""
import argparse
from pathlib import Path
import local_j4

def run(a):
    original=local_j4.environment
    def env(**kwargs):
        e=original(**kwargs)
        e.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',WAVETRISK_STAGE192_SPARSE_OUTER='0',
                 WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1')
        if a.observer:e.update(WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200='0',WAVETRISK_PMU201='0',WAVETRISK_LIFECYCLE205=str(a.mode))
        return e
    local_j4.environment=env
    a.profile=False
    try:local_j4.run(a)
    finally:local_j4.environment=original
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('build','fixture','out'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--binary',required=True);p.add_argument('--ranks',type=int,default=4);p.add_argument('--oracles',action='store_true')
    p.add_argument('--observer',action='store_true')
    p.add_argument('--mode',type=int,choices=(0,1),default=1)
    run(p.parse_args())
