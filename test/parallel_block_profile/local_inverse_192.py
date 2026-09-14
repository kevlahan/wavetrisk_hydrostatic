#!/usr/bin/env python3
"""Run a local J4 Stage 192 case with recorded, explicit runtime selectors."""
import argparse
from pathlib import Path
import local_j4


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--binary',required=True)
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--ranks',type=int,default=4)
    p.add_argument('--control',action='store_true')
    p.add_argument('--oracles',action='store_true')
    p.add_argument('--profile',action='store_true')
    args=p.parse_args()
    original=local_j4.environment
    def configured(**kwargs):
        env=original(**kwargs)
        env.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',
                   WAVETRISK_STAGE192_SPARSE_OUTER='0' if args.control else '1')
        return env
    local_j4.environment=configured
    try:
        local_j4.run(args)
    finally:
        local_j4.environment=original


if __name__=='__main__':main()
