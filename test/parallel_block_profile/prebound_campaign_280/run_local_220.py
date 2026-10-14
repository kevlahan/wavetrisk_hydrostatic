"""Full adaptive/restart correctness with the accepted runtime selectors."""
import argparse
from pathlib import Path
from types import SimpleNamespace
import local_j4
SELECTORS=dict(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',
 WAVETRISK_STAGE192_SPARSE_OUTER='0',WAVETRISK_STAGE193_SCALAR_COLUMNS='1',
 WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1')
original_environment=local_j4.environment
def environment(**kw):
 env=original_environment(**kw);env.update(SELECTORS);return env
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('build','fixture','out'):p.add_argument('--'+name,type=Path,required=True)
 p.add_argument('--kind',choices=['candidate','checked'],required=True);p.add_argument('--oracles',action='store_true')
 a=p.parse_args();local_j4.environment=environment
 local_j4.run(SimpleNamespace(build=a.build,fixture=a.fixture,out=a.out,ranks=4,
 binary='bin-220-'+a.kind+'/climate',oracles=a.oracles,profile=False))
