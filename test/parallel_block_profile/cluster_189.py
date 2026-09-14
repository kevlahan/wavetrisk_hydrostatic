#!/usr/bin/env python3
"""Stage 189: exact Stage 188/candidate and checked/optimized cluster gates.

All supplied binaries must use FP_CONTRACT=off; retain build identities in
BUILD_NOTES. This wrapper reuses the established 83-task timing campaign.
"""
import argparse
import json
from pathlib import Path

from cluster_187 import exact_gate, run as run_campaign, save


def finish_cross_mode(out):
    path = out / 'results.json'
    result = json.loads(path.read_text())
    if result['status'] != 'passed':
        raise ValueError('Underlying restart/timing campaign did not pass')
    result['status'] = 'checking_cross_mode'
    result['comparison_mode'] = 'Stage 188 versus Stage 189; all builds contraction off'
    result['limitations'][0] = 'Checked/optimized equality tested on this fixture; no new cluster legacy comparison.'
    save(path, result)
    try:
        for label in ('baseline-warmup', 'candidate-warmup'):
            gate = exact_gate(out / 'checked-off', out / label, result['parameters']['run_id'], [4])
            gate['reference'] = str(out / 'checked-off')
            result['gates']['cross-mode-' + label] = gate
            if not gate['passed']:
                raise ValueError('Checked/optimized checkpoint mismatch: ' + label)
        result['status'] = 'passed'
        save(path, result)
    except BaseException as error:
        result.update(status='failed', error=repr(error))
        save(path, result)
        raise
    print('PASS Stage 189: restart, oracle, baseline/candidate and checked/optimized exact gates', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'candidate', 'checked', 'fixture', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--build-notes', required=True)
    parser.add_argument('--profiles', action='store_true')
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.ranks = 83
    args.repeats = 3
    args.expected_checkpoints = [4]
    args.require_post_restart_remap = True
    args.arithmetic_control = False
    run_campaign(args)
    finish_cross_mode(args.out)


if __name__ == '__main__':
    main()
