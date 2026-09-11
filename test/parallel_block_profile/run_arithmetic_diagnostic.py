#!/usr/bin/env python3
"""Compare J5 solver arithmetic controls on the original two-restart fixture.

Supply separately built executables; this driver never changes compiler flags,
sources or numerical tolerances. Runs are sequential diagnostic observations,
not performance measurements. The checked reference must use identical input.
"""
import argparse
import json
from pathlib import Path
import re

from checkpoint_gate import compare_runs
from profile_protocol import execute, fixture_identity, parameters, replace_parameter


def coverage(root):
    log = (root / 'run.log').read_text()
    events = [(int(m[1]) if m[1] else None) for m in re.finditer(
        r'Restarting from checkpoint\s+(\d+)|Remapping vertical coordinates', log)]
    if [e for e in events if e is not None] != [4, 5, 6]:
        raise ValueError('Expected checkpoint reloads 4, 5, 6: ' + str(root))
    for cp in (5, 6):
        start = events.index(cp) + 1
        end = next((i for i in range(start, len(events)) if events[i] is not None), len(events))
        if None not in events[start:end]:
            raise ValueError('Missing remap following checkpoint ' + str(cp))
    return {'reloads': [4, 5, 6], 'remap_after_each_new_restart': True,
            'steps': len(re.findall(r'^\d+\s+\S+\s+d dt', log, re.M))}


def run(args):
    out = args.out.resolve()
    fixture = args.fixture.resolve(strict=True)
    reference = args.checked_reference.resolve(strict=True)
    frozen = fixture_identity(fixture)
    if frozen['checkpoint'] != 'test_checkpoint_0004.bin.zst':
        raise ValueError('Expected original J5 checkpoint-4 fixture')
    reference_identity = json.loads((reference / 'identity.json').read_text())
    if reference_identity['ranks'] != args.ranks:
        raise ValueError('Checked reference uses a different MPI rank count')
    expected = parameters(replace_parameter((fixture / 'simple.in').read_text(),
                                           'time_end', args.time_end))
    if parameters((reference / 'simple.in').read_text()) != expected or fixture_identity(reference) != frozen:
        raise ValueError('Checked reference uses different input, seed or grids')
    out.mkdir(parents=True, exist_ok=False)
    result = {'status': 'running', 'fixture': frozen, 'runs': {}, 'gates': {},
              'checked_reference_identity': reference_identity,
              'scope': 'Exact J5 checkpoint comparisons; no timing or universal equivalence claim.'}
    def save():
        (out / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    try:
        save()
        for label, binary in [('block', args.block), ('legacy', args.legacy)]:
            if fixture_identity(fixture) != frozen:
                raise ValueError('Fixture changed during diagnostic')
            result['active'] = label
            save()
            record = execute(binary.resolve(strict=True), fixture, out / label,
                             args.ranks, args.time_end, oracles=False, monitor=True)
            record['coverage'] = coverage(out / label)
            result['runs'][label] = record
            save()
        for root in (reference, out / 'block', out / 'legacy'):
            if parameters((root / 'simple.in').read_text()) != parameters((reference / 'simple.in').read_text()):
                raise ValueError('Unmatched input parameters: ' + str(root))
            if fixture_identity(root) != frozen:
                raise ValueError('Unmatched seed/grid identity: ' + str(root))
            coverage(root)
        result['gates']['checked_vs_block'] = compare_runs(
            reference, out / 'block', frozen['checkpoint'], exact=True)
        result['gates']['legacy_vs_block'] = compare_runs(
            out / 'legacy', out / 'block', frozen['checkpoint'], exact=True)
        result['gates']['legacy_screen'] = compare_runs(
            out / 'legacy', out / 'block', frozen['checkpoint'], exact=False)
        result['status'] = ('passed' if all(g['passed'] for g in result['gates'].values())
                            else 'completed_with_differences')
        result.pop('active', None)
        save()
        print(json.dumps({k: v['passed'] for k, v in result['gates'].items()}), flush=True)
    except Exception as error:
        result.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('block', 'legacy', 'checked-reference', 'fixture', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--time-end', default='0.3350')
    run(parser.parse_args())
