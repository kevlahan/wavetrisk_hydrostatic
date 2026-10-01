"""Independently verify the sealed Stage251 feature-ladder return."""
import hashlib
import json
from pathlib import Path
import statistics
import sys
import tarfile


STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
CANDIDATE_CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
MAIN_CHECKPOINT = '5a60b279f97aa48098126201e4e5a549f759a04c8a9d0ef18e94648ccba49765'
EARLY = ('BLOCK_COLUMNS', 'BLOCK_PRIMITIVES', 'BLOCK_DIFFUSION',
         'BLOCK_MASS_OWNER', 'BLOCK_VELOCITY_OWNER', 'BLOCK_TENDENCY_INITIAL')
LATE = ('BLOCK_RK', 'PHYSICS_COLUMN_LOCAL', 'PHYSICS_VELOCITY_PLAN', 'PHYSICS_OWNER',
        'PHYSICS_OWNER_BATCH', 'PHYSICS_NATIVE_INPUT', 'PHYSICS_NATIVE_STATE',
        'SOIL_NATIVE', 'ADAPT_NATIVE', 'LIFECYCLE_NATIVE', 'COEFFICIENT_NATIVE',
        'MIRROR_FREE')


def equal(actual, expected):
    if abs(actual-expected) > 1e-11:
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
        if result['status'] != 'complete' or result['decision'] != 'review_feature_ladder':
            raise ValueError('Campaign did not complete')
        if result['host'] != 'bb13.rhpcs.internal' or result['slurm_job_id'] != '55256':
            raise ValueError('Host or job differs')
        builds = result['build']['builds']
        if set(builds) != {'main'} | {f'c{i}' for i in range(7)}:
            raise ValueError('Build ladder differs')
        if len(builds['main']['physics_objects_sha256']) != 23:
            raise ValueError('Shared physics objects missing')
        for stage in range(7):
            command = builds[f'c{stage}']['command']
            selectors = dict(item.split('=', 1) for item in command if '=' in item)
            if command[:2] != ['make', '-j1'] or selectors['FP_CONTRACT'] != 'off':
                raise ValueError('Candidate build method differs')
            if selectors['STAGE202'] != 'on' or selectors['STAGE203'] != 'on':
                raise ValueError('Common block foundation differs')
            if any(selectors[name] != ('on' if i < stage else 'off') for i, name in enumerate(EARLY)):
                raise ValueError('Early-feature ladder differs')
            if any(selectors[name] != 'off' for name in LATE):
                raise ValueError('Later-feature selector differs')
        if 'ffp-contract=off' not in stream.extractfile(members['results/build-main.log']).read().decode():
            raise ValueError('Main FP setting differs')
        plan = read_json('placement-83.json')
        probe = read_json('placement-83-verified.json')
        mapping = {item['rank']: item['cpu'] for item in plan['mapping']}
        if set(mapping) != set(range(83)) or len(set(mapping.values())) != 83:
            raise ValueError('Core placement differs')
        if plan['socket_rank_counts'] != [42, 41] or {item['rank'] for item in probe} != set(mapping):
            raise ValueError('Placement probe differs')
        for item in probe:
            if item['hostname'] != result['host'] or item['affinity'] != [mapping[item['rank']]]:
                raise ValueError('Probe affinity differs')

        schedule = result['schedule']
        expected = ([f'J5-n83-validation-{kind}-{kind}' for kind in ['main']+[f'c{i}' for i in range(7)]]
                    + [f'J5-n83-confirm-{batch}-{index}-{kind}'
                       for batch in (1, 2)
                       for index, kind in enumerate(('c4', 'c5', 'c5', 'c4'), 1)])
        if schedule != expected or set(result['runs']) != set(schedule) or len(result['gates']) != 15:
            raise ValueError('Run schedule or gate count differs')
        medians = {}
        for name in schedule:
            row = result['runs'][name]
            kind = row['kind']
            if (row['returncode'] or row['ranks'] != 83 or row['steps'] != 12 or
                    row['remaps'] != 2 or row['vtk_members'] != 62 or
                    row['binary_sha256'] != builds[kind]['sha256'] or
                    row['state_sha256'] != STATE or
                    row['checkpoint_sha256'] != (MAIN_CHECKPOINT if kind == 'main' else CANDIDATE_CHECKPOINT)):
                raise ValueError('Run integrity differs: '+name)
            ordinary = row['ordinary']
            if ([item['index'] for item in ordinary['included']] != [1,2,3,4,5,6,9,10,11]
                    or [item['index'] for item in ordinary['excluded']] != [0,7,8]):
                raise ValueError('Ordinary-step selection differs')
            median = statistics.median(item['seconds'] for item in ordinary['included'])
            equal(median, ordinary['median'])
            medians[name] = median
            for rank in range(83):
                record = read_json(f'{name}/rank-{rank}.json')
                if record['returncode'] or record['hostname'] != result['host'] or record['affinity'] != [mapping[rank]]:
                    raise ValueError('Rank exit or binding differs: '+name)
        for name, gate in result['gates'].items():
            report = gate['report']
            if name not in schedule or result['runs'][name]['kind'] == 'main' or not gate['passed'] or gate['failures']:
                raise ValueError('Numerical gate differs: '+name)
            if any(report[key] for key in ('header_differences', 'topology_differences',
                                           'threshold_max_abs', 'threshold_nonfinite')):
                raise ValueError('Checkpoint threshold differs: '+name)

        validation = {kind: medians[f'J5-n83-validation-{kind}-{kind}']
                      for kind in ['main']+[f'c{i}' for i in range(7)]}
        recorded = result['decisions']['ladder']
        for kind, value in validation.items():
            equal(value, recorded['validation_medians'][kind])
        order = list(validation)
        adjacent = [{ 'before': order[i-1], 'after': order[i],
                      'difference_seconds': validation[order[i]]-validation[order[i-1]]}
                    for i in range(1, len(order))]
        for actual, claimed in zip(adjacent, recorded['adjacent']):
            if (actual['before'], actual['after']) != (claimed['before'], claimed['after']):
                raise ValueError('Adjacent stages differ')
            equal(actual['difference_seconds'], claimed['difference_seconds'])
        selected = max(adjacent, key=lambda row: row['difference_seconds'])
        if (selected['before'], selected['after']) != ('c4', 'c5'):
            raise ValueError('Largest jump differs')
        equal(selected['difference_seconds'], recorded['selected_adjacent']['difference_seconds'])
        before = [medians[name] for name in schedule[8:] if result['runs'][name]['kind'] == 'c4']
        after = [medians[name] for name in schedule[8:] if result['runs'][name]['kind'] == 'c5']
        if len(before) != 4 or len(after) != 4:
            raise ValueError('Confirmation run count differs')
        confirmation = result['decisions']['confirmation']
        equal(statistics.mean(before), confirmation['before_mean_seconds'])
        equal(statistics.mean(after), confirmation['after_mean_seconds'])
        equal(statistics.mean(after)-statistics.mean(before), confirmation['difference_seconds'])
        pairs = [medians[schedule[i+1]]-medians[schedule[i]] for i in (8,12)]
        pairs += [medians[schedule[i+2]]-medians[schedule[i+3]] for i in (8,12)]
        pairs = [pairs[0], pairs[2], pairs[1], pairs[3]]
        for actual, claimed in zip(pairs, confirmation['paired_differences_seconds']):
            equal(actual, claimed)

    with path.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    return dict(status='verified', job_id=result['slurm_job_id'],
                archive_sha256=digest, verified_files=len(manifest),
                validation_seconds=validation, adjacent_seconds=adjacent,
                confirmation=dict(before='c4', after='c5', before_mean_seconds=statistics.mean(before),
                                  after_mean_seconds=statistics.mean(after),
                                  difference_seconds=statistics.mean(after)-statistics.mean(before),
                                  paired_differences_seconds=pairs,
                                  independent_runs_per_variant=4),
                numerical_gates=len(result['gates']))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_251.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]), indent=2))
