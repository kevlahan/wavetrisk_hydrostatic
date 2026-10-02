"""Independently review the sealed Stage255 owner depth-index trial return."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import tarfile

from analyze import STEP

STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
BINARIES = {'baseline': 'e3feffb3eaddfbf49c2745cac84807f02eceabecba42eb2967133781e258fa91',
            'optimized': '8ae771aa5074ff25181eff68b6febbea3963bb67e8e09e53c951a07683686b98'}
OBJECT_MAP_SHA = '9717172ac48acbe10744240b539d1f562f65bc4d230a6fc2228f64945857f031'
PHASES = ('complete timestep', 'native velocity source', 'native velocity gradient')
PROFILE_ROW = re.compile(r'^\s*(.+?):\s*(\d+)\s*\|\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*\|')


def equal(actual, expected):
    if abs(actual-expected) > 1e-10:
        raise ValueError(f'Recomputed {actual} differs from recorded {expected}')


def review(path):
    path = Path(path).resolve(strict=True)
    with tarfile.open(path) as stream:
        listed = [item for item in stream if item.isfile()]
        members = {item.name: item for item in listed}
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

        def read_text(name):
            return stream.extractfile(members['results/'+name]).read().decode()

        result = read_json('results.json')
        if (result['status'] != 'complete' or result['decision'] != 'review_before_promotion' or
                result['job_id'] != '55261' or result['stage251_job_id'] != '55256' or
                result['host'] != 'bb13.rhpcs.internal' or
                result['prior_results_sha256'] != '138ec179faea1cea0ed4c8afb1f6b126f5982b5d64e49cf39de8c67014565b51'):
            raise ValueError('Campaign identity or completion differs')
        generated = result['regenerated_physics_files']
        if (len(generated) != 23 or len(set(generated)) != 23 or
                any(not name.startswith('source/src/physics/simple_physics/phyparam/tmp/') for name in generated)):
            raise ValueError('Unexpected generated physics files')
        build = result['build']
        if (build['source_sha256'] != 'a65f465f483e789f68f0f454bb91088fe81309d4cf51d62950a0d5bdd299f2dc' or
                build['baseline_binary_sha256'] != BINARIES['baseline'] or
                build['optimized_binary_sha256'] != BINARIES['optimized'] or
                len(build['physics_objects']) != 23):
            raise ValueError('Build identity differs')
        objects_sha = hashlib.sha256(json.dumps(build['physics_objects'], sort_keys=True,
                                               separators=(',', ':')).encode()).hexdigest()
        if objects_sha != OBJECT_MAP_SHA:
            raise ValueError('Compiled physics objects differ')
        command = build['command']
        if (command[:2] != ['make', '-j1'] or 'PARAM=param_J5' not in command or
                'BLOCK_VELOCITY_OWNER=on' not in command or 'BLOCK_TENDENCY_INITIAL=off' not in command or
                'BUILD_DIR=build-rk255-opt-c5' not in command or 'BIN_DIR=bin-rk255-opt-c5' not in command):
            raise ValueError('Optimized c5 build selectors differ')
        if 'COMPLETE: depth-index optimization trial' not in read_text('slurm-job.log'):
            raise ValueError('Slurm completion log differs')
        placement = read_json('placement-83.json')
        mapping = {item['rank']: item['cpu'] for item in placement['mapping']}
        if (set(mapping) != set(range(83)) or len(set(mapping.values())) != 83 or
                placement['socket_rank_counts'] != [42, 41]):
            raise ValueError('Physical-core map differs')
        expected = [f'J5-n83-{mode}-{index}-{kind}' for mode in ('off', 'profile')
                    for index, kind in enumerate(('baseline', 'optimized', 'optimized', 'baseline'), 1)]
        if result['schedule'] != expected or set(result['runs']) != set(expected):
            raise ValueError('A-B-B-A schedule differs')

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
            analysis = read_json(name+'/analysis.json')
            log = read_text(name+'/run.log')
            if hashlib.sha256(log.encode()).hexdigest() != row['log_sha256']:
                raise ValueError('Run log hash differs: '+name)
            states = [re.sub(r'\s+cpu =.*', '', line.strip()) for line in log.splitlines() if STEP.match(line)]
            if states != [item['state'] for item in analysis['step_records']] or len(states) != 12:
                raise ValueError('Raw step states differ: '+name)
            if hashlib.sha256('\n'.join(states).encode()).hexdigest() != STATE:
                raise ValueError('Raw state hash differs: '+name)
            if log.count('Remapping vertical coordinates') != 2:
                raise ValueError('Raw remap count differs: '+name)
            ordinary = row['ordinary']
            if (ordinary != analysis['ordinary'] or
                    [item['index'] for item in ordinary['included']] != [1,2,3,4,5,6,9,10,11] or
                    [item['index'] for item in ordinary['excluded']] != [0,7,8]):
                raise ValueError('Ordinary-step window differs: '+name)
            medians[name] = statistics.median(item['seconds'] for item in ordinary['included'])
            equal(medians[name], ordinary['median'])
            for rank in range(83):
                record = read_json(f'{name}/rank-{rank}.json')
                if (record['rank'] != rank or record['returncode'] or
                        record['hostname'] != result['host'] or record['affinity'] != [mapping[rank]]):
                    raise ValueError('Rank placement/exit differs: '+name)
                if f'results/{name}/rank-{rank}.resources' not in members:
                    raise ValueError('Rank resource file absent: '+name)
            if mode == 'off':
                if 'Parallel-block profiling summary (inclusive):' in log or row['profile'] is not None:
                    raise ValueError('Unprofiled run emitted a profile: '+name)
                continue
            if log.count('Parallel-block profiling summary (inclusive):') != 2:
                raise ValueError('Profile-window count differs: '+name)
            found = {phase: [] for phase in PHASES}
            for line in log.splitlines():
                match = PROFILE_ROW.match(line)
                if match and match[1] in found:
                    found[match[1]].append((int(match[2]), float(match[4])))
            if any(len(values) != 2 for values in found.values()):
                raise ValueError('Profile phase coverage differs: '+name)
            phase_totals[name] = {}
            for phase, values in found.items():
                recorded = row['profile'][phase]
                if [item['calls'] for item in recorded['windows']] != [value[0] for value in values]:
                    raise ValueError('Profile call count differs: '+name)
                for index, item in enumerate(recorded['windows']):
                    equal(item['seconds_avg'], values[index][1])
                total = sum(value[1] for value in values)
                equal(total, recorded['seconds_avg_sum'])
                phase_totals[name][phase] = total
            if ([v[0] for v in found['complete timestep']] != [8,4] or
                    [v[0] for v in found['native velocity source']] != [2880,1320] or
                    [v[0] for v in found['native velocity gradient']] != [960,480]):
                raise ValueError('Profile call counts differ: '+name)

        decisions = {}
        for mode in ('off', 'profile'):
            values = [medians[f'J5-n83-{mode}-{i}-{kind}']
                      for i, kind in enumerate(('baseline', 'optimized', 'optimized', 'baseline'), 1)]
            baseline = statistics.mean((values[0], values[3]))
            optimized = statistics.mean((values[1], values[2]))
            recorded = result['decisions'][mode]
            for actual, key in ((baseline, 'baseline_mean_seconds'),
                                (optimized, 'optimized_mean_seconds'),
                                (optimized-baseline, 'difference_seconds')):
                equal(actual, recorded[key])
            decisions[mode] = dict(baseline_mean_seconds=baseline,
                                   optimized_mean_seconds=optimized,
                                   difference_seconds=optimized-baseline,
                                   paired_differences_seconds=[values[1]-values[0], values[2]-values[3]])
        phase_deltas = {}
        for phase in PHASES:
            baseline = statistics.mean(phase_totals[f'J5-n83-profile-{i}-baseline'][phase] for i in (1,4))
            optimized = statistics.mean(phase_totals[f'J5-n83-profile-{i}-optimized'][phase] for i in (2,3))
            phase_deltas[phase] = dict(baseline_total_seconds=baseline, optimized_total_seconds=optimized,
                                       difference_seconds=optimized-baseline,
                                       difference_per_step_seconds=(optimized-baseline)/12)

    with path.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    return dict(status='verified', job_id=result['job_id'], archive_sha256=digest,
                verified_files=len(manifest), full_j5_runs=8,
                generated_physics_files=len(generated), decision=decisions,
                phase_deltas=phase_deltas)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_255.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]), indent=2))
