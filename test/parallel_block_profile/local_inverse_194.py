#!/usr/bin/env python3
"""Run a local J4 Stage 194 case with recorded, explicit runtime selectors."""
import argparse
import json
import re
from analyze_scalar_inverse_193 import scalar_work
from analyze_inner_inverse_194 import inner_work, split_profile
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
                   WAVETRISK_STAGE192_SPARSE_OUTER='0',
                   WAVETRISK_STAGE193_SCALAR_COLUMNS='1',
                   WAVETRISK_STAGE194_INNER_GEOMETRY='0' if args.control else '1')
        return env
    local_j4.environment=configured
    try:
        local_j4.run(args)
        log=(args.out/'run.log').read_text()
        selected='0' if args.control else '1'
        if re.findall(r'Stage 194 requested inner inverse geometry = ([01])',log)!=[selected]:
            raise ValueError('Missing or incorrect Stage 194 selector: reject stale binary')
        if re.findall(r'Stage 192 requested sparse outer install = ([01])',log)!=['0']:
            raise ValueError('Stage 192 sparse transfer must remain disabled')
        if re.findall(r'Stage 193 requested scalar inverse columns = ([01])',log)!=['1']:
            raise ValueError('Validated Stage 193 scalar kernel must remain enabled')
        if args.profile:
            report={'work':inner_work(args.out/'run.log'), 'timing':split_profile(args.out/'run.log')}
            (args.out/'inner-inverse.json').write_text(json.dumps(report,indent=2)+'\n')
        if args.profile and not args.control:
            (args.out/'scalar-work.json').write_text(json.dumps(scalar_work(args.out/'run.log'),indent=2)+'\n')
    except BaseException as error:
        path=args.out/'results.json'
        if path.exists():
            result=json.loads(path.read_text())
            result.update(status='failed',error=repr(error))
            path.write_text(json.dumps(result,indent=2)+'\n')
        raise
    finally:
        local_j4.environment=original


if __name__=='__main__':main()
