"""Stage222 cleanup correctness: use accepted selectors and adaptive/restart fixture."""
import argparse
from pathlib import Path
from types import SimpleNamespace
import local_j4
from run_local_220 import SELECTORS

original_environment = local_j4.environment


def environment(**kw):
    env = original_environment(**kw)
    env.update(SELECTORS)
    return env


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'fixture', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--binary', required=True)
    p.add_argument('--oracles', action='store_true')
    a = p.parse_args()
    local_j4.environment = environment
    local_j4.run(SimpleNamespace(build=a.build, fixture=a.fixture, out=a.out,
                                ranks=4, binary=a.binary, oracles=a.oracles, profile=False))
