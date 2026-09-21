"""Local checked prototype with a deliberately empty MPI rank zero."""
import argparse
from pathlib import Path
from types import SimpleNamespace
import local_j4
from run_local_220 import SELECTORS

original = local_j4.environment


def environment(**kw):
    env = original(**kw)
    env.update(SELECTORS, WAVETRISK_STAGE223_ALLOW='1', WAVETRISK_STAGE223_LAYOUT='empty-root')
    return env


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('build','fixture','out'):
        p.add_argument('--'+key,type=Path,required=True)
    a = p.parse_args()
    local_j4.environment = environment
    local_j4.run(SimpleNamespace(build=a.build,fixture=a.fixture,out=a.out,ranks=4,
                                binary='bin-223-checked/climate',oracles=False,profile=False))
