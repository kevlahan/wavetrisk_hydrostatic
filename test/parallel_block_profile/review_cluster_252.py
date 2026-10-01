"""Independently review the sealed Stage252 velocity phase-profile return."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import tarfile


STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
BINARIES = {'c4': '02ce866f8c82763eb5b600ea5563e7ce484a2e7f88fb63a1f3b8f161778712ec',
            'c5': 'e3feffb3eaddfbf49c2745cac84807f02eceabecba42eb2967133781e258fa91'}
PHASES = ('complete timestep', 'native velocity source', 'native velocity gradient')
PROFILE_ROW = re.compile(r'^\s*(.+?):\s*(\d+)\s*\|\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*\|')


def equal(actual, expected):
    if abs(actual-expected) > 1e-10:
        raise ValueError(f'Recomputed {actual} differs from recorded {expected}')


def review(path):
    path = Path(path).resolve(strict=True)
    with tarfile.open(path) as stream:
        listed = [member for member in stream if member.isfile()]
        members = {member.name: member for member in listed}
        if len(members) != len(listed):
            raise ValueError('Duplicate archive member')
        manifest = json.load(stream.extractfile(members['results/REVIEW_SHA256.json']))
        if set(members) != {'results/'+name for name in manifest} | {'results/REVIEW_SHA256.json'}:
            raise ValueError('Review inventory differs')
        for name, expected in manifest.items():
            if hashlib.sha256(stream.extractfile(members['results/'+name]).read()).hexdigest() != expected:
                raise ValueError('Review hash differs: '+name)

        def read_json(name):
            return json.load(stream.extractfile(members['results/'+name]))

        result = read_json('results.json')
        if (result['status'] != 'complete' or result['decision'] != 'review_profile_attribution' or
                result['job_id'] != '55258' or result['stage251_job_id'] != '55256' or
                result['host'] != 'bb13.rhpcs.internal'):
            raise ValueError('Campaign identity or completion differs')
        if result['prior_results_sha256'] != '138ec179faea1cea0ed4c8afb1f6b126f5982b5d64e49cf39de8c67014565b51':
            raise ValueError('Prior Stage251 result differs')
        plan = read_json('placement-83.json')
        mapping = {item['rank']: item['cpu'] for item in plan['mapping']}
        if set(mapping) != set(range(83)) or len(set(mapping.values())) != 83 or plan['socket_rank_counts'] != [42, 41]:
            raise ValueError('Core map differs')
        expected = [f'J5-n83-{mode}-{index}-{kind}' for mode in ('off', 'profile')
                    for index, kind in enumerate(('c4', 'c5', 'c5', 'c4'), 1)]
        if result['schedule'] != expected or set(result['runs']) != set(expected):
            raise ValueError('Run schedule differs')

        medians = {}
        phase_totals = {}
        for name in expected:
            row = result['runs'][name]
            kind, mode = row['kind'], row['mode']
            if (row['returncode'] or row['ranks'] != 83 or row['steps'] != 12 or
                    row['remaps'] != 2 or row['vtk_members'] != 62 or
                    row['binary_sha256'] != BINARIES[kind] or row['state_sha256'] != STATE or
                    row['checkpoint_sha256'] != CHECKPOINT):
                raise ValueError('Run science/binary gate differs: '+name)
            ordinary = row['ordinary']
            if ([item['index'] for item in ordinary['included']] != [1,2,3,4,5,6,9,10,11] or
                    [item['index'] for item in ordinary['excluded']] != [0,7,8]):
                raise ValueError('Ordinary-step window differs: '+name)
            medians[name] = statistics.median(item['seconds'] for item in ordinary['included'])
            equal(medians[name], ordinary['median'])
            for rank in range(83):
                record = read_json(f'{name}/rank-{rank}.json')
                if (record['rank'] != rank or record['returncode'] or
                        record['hostname'] != result['host'] or record['affinity'] != [mapping[rank]]):
                    raise ValueError('Rank placement/exit differs: '+name)
            log = stream.extractfile(members[f'results/{name}/run.log']).read().decode()
            if mode == 'off':
                if 'Parallel-block profiling summary (inclusive):' in log or row['profile'] is not None:
                    raise ValueError('Unprofiled run emitted profile: '+name)
                continue
            if log.count('Parallel-block profiling summary (inclusive):') != 2:
                raise ValueError('Wrong number of profile windows: '+name)
            found = {phase: [] for phase in PHASES}
            for line in log.splitlines():
                match = PROFILE_ROW.match(line)
                if match and match[1] in found:
                    found[match[1]].append((int(match[2]), float(match[3]),
                                            float(match[4]), float(match[5])))
            if any(len(values) != 2 for values in found.values()):
                raise ValueError('Incomplete profile phases: '+name)
            phase_totals[name] = {}
            for phase, values in found.items():
                recorded = row['profile'][phase]
                if [item['calls'] for item in recorded['windows']] != [v[0] for v in values]:
                    raise ValueError('Profile call count differs: '+name)
                for i, item in enumerate(recorded['windows']):
                    for actual, key in zip(values[i][1:], ('seconds_min', 'seconds_avg', 'seconds_max')):
                        equal(actual, item[key])
                avg = sum(value[2] for value in values)
                equal(avg, recorded['seconds_avg_sum'])
                phase_totals[name][phase] = avg
            if ([v[0] for v in found['complete timestep']] != [8,4] or
                    [v[0] for v in found['native velocity source']] != [2880,1320] or
                    [v[0] for v in found['native velocity gradient']] != [960,480]):
                raise ValueError('Profile coverage differs: '+name)

        decisions = {}
        for mode in ('off', 'profile'):
            runs = [medians[f'J5-n83-{mode}-{i}-{kind}']
                    for i, kind in enumerate(('c4', 'c5', 'c5', 'c4'), 1)]
            c4 = statistics.mean((runs[0], runs[3]))
            c5 = statistics.mean((runs[1], runs[2]))
            recorded = result['decisions'][mode]
            for actual, key in ((c4, 'c4_mean_seconds'), (c5, 'c5_mean_seconds'),
                                (c5-c4, 'difference_seconds')):
                equal(actual, recorded[key])
            decisions[mode] = dict(c4_mean_seconds=c4, c5_mean_seconds=c5,
                                   difference_seconds=c5-c4,
                                   paired_differences_seconds=[runs[1]-runs[0], runs[2]-runs[3]])
        phase_deltas = {}
        for phase in PHASES:
            values = {kind: statistics.mean(phase_totals[f'J5-n83-profile-{i}-{kind}'][phase]
                                            for i in ((1,4) if kind == 'c4' else (2,3)))
                      for kind in ('c4', 'c5')}
            phase_deltas[phase] = dict(c4_total_seconds=values['c4'], c5_total_seconds=values['c5'],
                                       difference_seconds=values['c5']-values['c4'],
                                       difference_per_step_seconds=(values['c5']-values['c4'])/12)
        source_fraction = (phase_deltas['native velocity source']['difference_seconds'] /
                           phase_deltas['complete timestep']['difference_seconds'])

    with path.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    return dict(status='verified', job_id=result['job_id'], archive_sha256=digest,
                verified_files=len(manifest), full_j5_runs=8,
                decision=decisions, phase_deltas=phase_deltas,
                velocity_source_fraction_of_profiled_timestep_gap=source_fraction)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_252.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]), indent=2))
