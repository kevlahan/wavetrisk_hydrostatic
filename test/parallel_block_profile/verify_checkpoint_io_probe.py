#!/usr/bin/env python3
"""Read saved synthetic probe files after MPI exits; no NumPy/MPI required."""
import argparse
import json
from pathlib import Path
import re
import struct

DOMAINS = 160
PAYLOAD = 8192
DATA_POS = 24 + 20 * DOMAINS
HEADER = (0x5741564554524953, 1, DOMAINS)
DIRECTORY = (struct.pack('<160i', *([1] * DOMAINS)) +
             struct.pack('<160q', *(DATA_POS + 16 * i for i in range(DOMAINS))) +
             struct.pack('<160q', *([16] * DOMAINS)))


def verify(path, ranks):
    data = path.read_bytes()
    header = struct.unpack('<3q', data[:24]) if len(data) >= 24 else None
    failures = []
    for rank in range(ranks):
        start = DATA_POS + rank * PAYLOAD
        block = data[start:start + PAYLOAD]
        expected = rank % 100 + 1
        wrong = [i for i, value in enumerate(block) if value != expected]
        if wrong or len(block) != PAYLOAD:
            first = wrong[0] if wrong else len(block)
            failures.append({'rank_block': rank, 'file_offset': start + first,
                             'mismatched_or_missing_bytes': len(wrong) + PAYLOAD - len(block),
                             'expected_byte': expected, 'observed_byte': block[first] if first < len(block) else None})
    return {'file': str(path), 'bytes': len(data), 'header': header,
            'bad_size': len(data) != DATA_POS + ranks * PAYLOAD,
            'bad_header': header != HEADER, 'bad_directory': data[24:DATA_POS] != DIRECTORY,
            'bad_payload_blocks': failures}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--ranks', type=int, required=True)
    args = parser.parse_args()
    if args.ranks < 2:
        parser.error('at least two ranks required')
    files = sorted(args.directory.glob('probe-*-*.bin'))
    if not files:
        parser.error('no probe files found')
    reports, groups = [], {}
    for path in files:
        match = re.fullmatch(r'probe-(\d+)-(\d+)\.bin', path.name)
        if not match:
            continue
        report = verify(path, args.ranks)
        reports.append(report)
        group = groups.setdefault(match[1], {'files': 0, 'bad_size': 0, 'bad_header': 0,
                                              'bad_directory': 0, 'bad_payload_blocks': 0})
        group['files'] += 1
        for key in ('bad_size', 'bad_header', 'bad_directory'):
            group[key] += int(report[key])
        group['bad_payload_blocks'] += len(report['bad_payload_blocks'])
        if any(report[k] for k in ('bad_size', 'bad_header', 'bad_directory', 'bad_payload_blocks')):
            print('PERSISTENT FAILURE', json.dumps(report), flush=True)
    if not reports:
        parser.error('no recognized probe filenames')
    result = {'directory': str(args.directory.resolve()), 'ranks': args.ranks, 'groups': groups, 'files': reports}
    (args.directory / 'postrun-verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Post-run verification:', json.dumps(groups, sort_keys=True), flush=True)
    # Diagnostic campaigns keep running even when an earlier mode fails.
    # Counts, not this process exit code, determine whether a mode passed.


if __name__ == '__main__':
    main()
