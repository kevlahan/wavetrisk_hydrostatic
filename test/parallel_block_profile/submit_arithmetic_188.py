#!/usr/bin/env python3
"""Validate the recorded builds, then run the 83-task arithmetic campaign."""
import argparse
import json
from pathlib import Path
import shutil
import sys

from cluster_187 import digest, run


def main():
    if len(sys.argv) != 5:
        raise SystemExit('Usage: submit_arithmetic_188.py KIT TREE FIXTURE NEW_OUTPUT')
    kit, tree, fixture, out = [Path(p).resolve() for p in sys.argv[1:]]
    if out.exists():
        raise ValueError('Use a new output directory')
    identity = json.loads((tree / 'arithmetic-188-builds.json').read_text())
    binaries = {}
    for kind, record in identity['builds'].items():
        path = tree / record['binary']
        if digest(path) != record['sha256']:
            raise ValueError('Executable changed since build: ' + str(path))
        binaries[kind] = path
    args = argparse.Namespace(**binaries, fixture=fixture, out=out, ranks=83,
                              repeats=3, profiles=False, expected_checkpoints=[4],
                              require_post_restart_remap=True, arithmetic_control=True,
                              build_notes=json.dumps(identity, sort_keys=True))
    run(args)
    shutil.copy2(tree / 'arithmetic-188-builds.json', out / 'builds.json')


if __name__ == '__main__':
    main()
