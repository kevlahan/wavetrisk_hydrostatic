"""Independently review the sealed Stage249 direct-main J5 return archive."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import tarfile


EXPECTED_STATES = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
EXPECTED_CANDIDATE_CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
MAIN_COMMIT = '6a245db7f1632b0205c2bebe95e5f0cd535c2582'
SCHEDULE = [(16, 'validation', 'main'), (16, 'validation', 'candidate'),
            (83, 'validation', 'main'), (83, 'validation', 'candidate'),
            (83, 'confirm', 'main'), (83, 'confirm', 'candidate'),
            (83, 'confirm', 'candidate'), (83, 'confirm', 'main')]


def close(actual, expected):
    if abs(actual - expected) > 1e-11:
        raise ValueError(f'Calculated value {actual} differs from recorded {expected}')


def review(path):
    path = Path(path).resolve(strict=True)
    with tarfile.open(path) as stream:
        listed = [entry for entry in stream if entry.isfile()]
        files = {entry.name: entry for entry in listed}
        if len(files) != len(listed):
            raise ValueError('Duplicate review member')
        manifest = json.load(stream.extractfile(files['results/REVIEW_SHA256.json']))
        if set(files) != {'results/' + name for name in manifest} | {'results/REVIEW_SHA256.json'}:
            raise ValueError('Review member inventory differs')
        for name, digest in manifest.items():
            data = stream.extractfile(files['results/' + name]).read()
            if hashlib.sha256(data).hexdigest() != digest:
                raise ValueError('Review member checksum differs: ' + name)

        def read_json(name):
            return json.load(stream.extractfile(files['results/' + name]))

        result = read_json('results.json')
        if result['status'] != 'complete' or result['decision'] != 'review_regression':
            raise ValueError('Direct comparison did not complete with expected disposition')
        if result['host'] != 'bb13.rhpcs.internal' or result['main_identity']['main_commit'] != MAIN_COMMIT:
            raise ValueError('Wrong node or main source')
        if result['build']['main_identity'] != result['main_identity']:
            raise ValueError('Build and run main identities differ')
        if len(result['build']['builds']['main']['physics_objects_sha256']) != 23:
            raise ValueError('Incomplete shared physics build')
        if len(result['schedule']) != len(SCHEDULE) or len(result['runs']) != len(SCHEDULE):
            raise ValueError('Incomplete direct comparison schedule')
        if [(result['runs'][name]['ranks'], result['runs'][name]['phase'],
             result['runs'][name]['kind']) for name in result['schedule']] != SCHEDULE:
            raise ValueError('Run order is not validation followed by A-B-B-A')

        placements = {}
        for ranks in (16, 83):
            planned = read_json(f'placement-{ranks}.json')
            probed = read_json(f'placement-{ranks}-verified.json')
            mapping = {row['rank']: row['cpu'] for row in planned['mapping']}
            if set(mapping) != set(range(ranks)) or len(set(mapping.values())) != ranks:
                raise ValueError('Placement has missing or repeated cores')
            if {row['rank'] for row in probed} != set(mapping):
                raise ValueError('Placement probe has missing ranks')
            for row in probed:
                if row['hostname'] != result['host'] or row['affinity'] != [mapping[row['rank']]]:
                    raise ValueError('Placement probe differs from plan')
            placements[str(ranks)] = planned['socket_rank_counts']

        runs = []
        checkpoints = {'main': set(), 'candidate': set()}
        for name in result['schedule']:
            row = result['runs'][name]
            ranks, kind = row['ranks'], row['kind']
            if row['returncode'] or row['steps'] != 12 or row['remaps'] != 2 or row['vtk_members'] != 62:
                raise ValueError('Incomplete run: ' + name)
            if row['state_sha256'] != EXPECTED_STATES:
                raise ValueError('Adaptive state differs: ' + name)
            checkpoints[kind].add(row['checkpoint_sha256'])
            if kind == 'candidate':
                if row['checkpoint_sha256'] != EXPECTED_CANDIDATE_CHECKPOINT:
                    raise ValueError('Candidate differs from prior qualified J5 result')
                if row['storage']['absence_checks'] != 14*ranks or row['storage']['bounded_output_checks'] != 2*ranks:
                    raise ValueError('Incomplete candidate storage evidence')
            log = stream.extractfile(files[f'results/{name}/run.log']).read().decode()
            absence = [int(x) for x in re.findall(r'Domain prognostic storage absent rank=(\d+)', log)]
            output = [int(x) for x in re.findall(r'Native output alias adapter rank=(\d+)', log)]
            staging = [int(x) for x in re.findall(r'Legacy scalar staging and RK snapshots absent rank=(\d+)', log)]
            if kind == 'candidate':
                if any(absence.count(rank) != 14 or output.count(rank) != 2 or staging.count(rank) != 2
                       for rank in range(ranks)):
                    raise ValueError('Candidate storage diagnostics differ: ' + name)
            elif absence or output or staging:
                raise ValueError('Main unexpectedly emitted candidate storage diagnostics')
            ordinary = row['ordinary']
            included = ordinary['included']
            excluded = ordinary['excluded']
            if [x['index'] for x in included] != [1, 2, 3, 4, 5, 6, 9, 10, 11]:
                raise ValueError('Wrong ordinary step window')
            if [x['index'] for x in excluded] != [0, 7, 8]:
                raise ValueError('Wrong excluded step window')
            median = statistics.median(x['seconds'] for x in included)
            close(median, ordinary['median'])
            rank_rows = [read_json(f'{name}/rank-{rank}.json') for rank in range(ranks)]
            mapping = read_json(f'placement-{ranks}.json')['mapping']
            for rank, rank_row in enumerate(rank_rows):
                if (rank_row['returncode'] or rank_row['hostname'] != result['host']
                        or rank_row['affinity'] != [mapping[rank]['cpu']]):
                    raise ValueError('Rank placement or exit differs: ' + name)
            runs.append(dict(name=name, kind=kind, ranks=ranks, phase=row['phase'],
                             median_seconds=median, checkpoint_sha256=row['checkpoint_sha256'],
                             max_rank_peak_rss_kib=row['max_rank_peak_rss_kib'],
                             sum_user_cpu_seconds=sum(x['user_cpu_seconds'] for x in rank_rows),
                             storage_diagnostic_lines=len(absence)+len(output)+len(staging)))
        if any(len(values) != 1 for values in checkpoints.values()):
            raise ValueError('Same-version checkpoint hashes differ')

        if len(result['gates']) != 4:
            raise ValueError('Missing cross-version numerical gates')
        for name, gate in result['gates'].items():
            report = gate['report']
            if not gate['passed'] or gate['failures'] or gate['exact_observed']:
                raise ValueError('Unexpected numerical disposition: ' + name)
            if any(report[field] for field in ('header_differences', 'topology_differences',
                                               'threshold_max_abs', 'threshold_nonfinite')):
                raise ValueError('Checkpoint header/topology/threshold differs: ' + name)
            for field, values in report['fields'].items():
                group = field.split('/')[1]
                medium = field.split('/')[2]
                bound = 0 if medium == 'soil' else gate['policy']['atmosphere_absolute'][group]
                if values['nonfinite'] or values['max_abs'] > bound:
                    raise ValueError('Checkpoint field exceeds numerical policy: ' + field)

        confirm = result['decisions']['confirmation']
        main = [runs[i]['median_seconds'] for i in (4, 7)]
        candidate = [runs[i]['median_seconds'] for i in (5, 6)]
        main_mean, candidate_mean = statistics.mean(main), statistics.mean(candidate)
        pair_fractions = [candidate[0]/main[0]-1, candidate[1]/main[1]-1]
        close(main_mean, confirm['baseline_run_median_mean'])
        close(candidate_mean, confirm['candidate_run_median_mean'])
        close(candidate_mean-main_mean, -confirm['gain_seconds'])
        close(candidate_mean/main_mean-1, confirm['regression_fraction'])
        for actual, expected in zip(pair_fractions, confirm['paired_regression_fractions']):
            close(actual, expected)
        if confirm['decision'] != 'review_regression' or confirm['independent_runs_per_version'] != 2:
            raise ValueError('Wrong direct comparison decision')

    with path.open('rb') as source:
        archive_sha = hashlib.file_digest(source, 'sha256').hexdigest()
    return dict(status='verified', job_id=result['slurm_job_id'], review_archive=str(path),
                review_archive_sha256=archive_sha, verified_files=len(manifest),
                main_commit=MAIN_COMMIT, host=result['host'], runs=runs, placements=placements,
                numerical_gates=len(result['gates']), candidate_matches_prior_j5=True,
                main_checkpoint_sha256=next(iter(checkpoints['main'])),
                candidate_checkpoint_sha256=next(iter(checkpoints['candidate'])),
                confirmation=dict(main_mean_seconds=main_mean, candidate_mean_seconds=candidate_mean,
                                  candidate_minus_main_seconds=candidate_mean-main_mean,
                                  candidate_to_main_ratio=candidate_mean/main_mean,
                                  candidate_regression_fraction=candidate_mean/main_mean-1,
                                  paired_regression_fractions=pair_fractions,
                                  independent_runs_per_version=2, decision=confirm['decision']))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_249.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]), indent=2))
