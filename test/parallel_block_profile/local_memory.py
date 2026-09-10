#!/usr/bin/env python3
"""Read-only macOS memory observations during an existing local experiment.

Counters are SYSTEM-WIDE; differences are not attributable to climate alone.
RSS is a sample, not peak RSS or compressed footprint. No simulation launches.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--experiment', type=Path, required=True)
parser.add_argument('--expected-runs', type=int, required=True)
parser.add_argument('--seconds', type=int, default=900)
args = parser.parse_args()
with (args.experiment / 'memory-observations.jsonl').open('x') as output:
    deadline = time.monotonic() + args.seconds
    while time.monotonic() < deadline:
        paths = sorted(p.name for p in args.experiment.iterdir() if p.is_dir())
        row = {'timestamp': time.time(), 'latest_run': paths[-1] if paths else None}
        for key, command in [('vm_stat', ['vm_stat']), ('swap', ['sysctl', 'vm.swapusage']),
                             ('processes', ['ps', '-Ao', 'pid,ppid,rss,pcpu,comm'])]:
            result = subprocess.run(command, capture_output=True, text=True)
            text = result.stdout
            if key == 'processes':
                text = '\n'.join(line for line in text.splitlines() if './climateJ5' in line)
            row[key] = text
            if result.returncode:
                row[key + '_error'] = result.stderr
        output.write(json.dumps(row) + '\n')
        output.flush()
        try:
            results = json.loads((args.experiment / 'results.json').read_text())
            if len(results['records']) >= args.expected_runs:
                break
        except (FileNotFoundError, json.JSONDecodeError):
            pass
        time.sleep(5)
