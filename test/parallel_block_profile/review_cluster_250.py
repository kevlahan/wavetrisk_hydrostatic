"""Independently review the sealed Stage250 logging-calibration return."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
import tarfile


STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
CANDIDATE_CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
MAIN_CHECKPOINT = '5a60b279f97aa48098126201e4e5a549f759a04c8a9d0ef18e94648ccba49765'
MODES = (['main', 'loud', 'quiet'] + ['loud', 'quiet', 'quiet', 'loud'] * 2
         + ['main', 'quiet', 'quiet', 'main'])
PHASES = (['validation'] * 3 + ['logging'] * 8 + ['residual'] * 4)


def equal(actual, expected):
    if abs(actual-expected) > 1e-11:
        raise ValueError(f'Recomputed {actual} differs from recorded {expected}')


def review(path):
    path = Path(path).resolve(strict=True)
    with tarfile.open(path) as stream:
        listed = [entry for entry in stream if entry.isfile()]
        members = {entry.name: entry for entry in listed}
        if len(members) != len(listed):
            raise ValueError('Duplicate archive member')
        manifest = json.load(stream.extractfile(members['results/REVIEW_SHA256.json']))
        if set(members) != {'results/'+name for name in manifest} | {'results/REVIEW_SHA256.json'}:
            raise ValueError('Review inventory differs')
        for name, expected in manifest.items():
            data = stream.extractfile(members['results/'+name]).read()
            if hashlib.sha256(data).hexdigest() != expected:
                raise ValueError('Review hash differs: '+name)

        def read_json(name):
            return json.load(stream.extractfile(members['results/'+name]))

        result = read_json('results.json')
        if result['status'] != 'complete' or result['decision'] != 'review_logging_and_residual':
            raise ValueError('Campaign did not complete')
        if result['host'] != 'bb13.rhpcs.internal' or len(result['schedule']) != 15:
            raise ValueError('Wrong host or run count')
        if len(result['build']['builds']['main']['physics_objects_sha256']) != 23:
            raise ValueError('Shared physics build incomplete')
        if result['build']['diagnostic_patch']['patched_domain_sha256'] != '637413d215228154c3c0987b2be4a8b1a6d815653aeaaf4236187fd53015e1f5':
            raise ValueError('Wrong diagnostic patch')
        plan = read_json('placement-83.json')
        probe = read_json('placement-83-verified.json')
        mapping = {row['rank']: row['cpu'] for row in plan['mapping']}
        if set(mapping) != set(range(83)) or len(set(mapping.values())) != 83:
            raise ValueError('Placement has missing or repeated cores')
        if plan['socket_rank_counts'] != [42, 41] or {row['rank'] for row in probe} != set(mapping):
            raise ValueError('Placement probe or socket counts differ')
        for row in probe:
            if row['hostname'] != result['host'] or row['affinity'] != [mapping[row['rank']]]:
                raise ValueError('Placement probe differs from plan')

        runs = []
        for name, mode, phase in zip(result['schedule'], MODES, PHASES):
            row = result['runs'][name]
            if row['mode'] != mode or row['phase'] != phase or row['ranks'] != 83:
                raise ValueError('Unexpected run schedule: '+name)
            if row['returncode'] or row['steps'] != 12 or row['remaps'] != 2 or row['vtk_members'] != 62:
                raise ValueError('Incomplete run: '+name)
            kind = 'main' if mode == 'main' else 'candidate'
            if row['kind'] != kind or row['binary_sha256'] != result['build']['builds'][kind]['sha256']:
                raise ValueError('Wrong binary: '+name)
            if row['state_sha256'] != STATE:
                raise ValueError('Adaptive trajectory differs: '+name)
            if row['checkpoint_sha256'] != (MAIN_CHECKPOINT if kind == 'main' else CANDIDATE_CHECKPOINT):
                raise ValueError('Checkpoint differs: '+name)
            ordinary = row['ordinary']
            included = ordinary['included']
            if [x['index'] for x in included] != [1,2,3,4,5,6,9,10,11]:
                raise ValueError('Ordinary timing window differs')
            if [x['index'] for x in ordinary['excluded']] != [0,7,8]:
                raise ValueError('Excluded timing window differs')
            median = statistics.median(x['seconds'] for x in included)
            equal(median, ordinary['median'])
            log = stream.extractfile(members[f'results/{name}/run.log']).read().decode()
            absence = [int(x) for x in re.findall(r'Domain prognostic storage absent rank=(\d+)',log)]
            output = [int(x) for x in re.findall(r'Native output alias adapter rank=(\d+)',log)]
            staging = [int(x) for x in re.findall(r'Legacy scalar staging and RK snapshots absent rank=(\d+)',log)]
            for rank in range(83):
                if (absence.count(rank) != (1162//83 if mode=='loud' else 0)
                        or output.count(rank) != (0 if mode=='main' else 2)
                        or staging.count(rank) != (0 if mode=='main' else 2)):
                    raise ValueError('Storage logging mode differs: '+name)
            for rank in range(83):
                rank_row = read_json(f'{name}/rank-{rank}.json')
                if (rank_row['returncode'] or rank_row['hostname'] != result['host']
                        or rank_row['affinity'] != [mapping[rank]]):
                    raise ValueError('Rank exit or placement differs: '+name)
            runs.append(dict(name=name, mode=mode, phase=phase, median_seconds=median,
                             diagnostic_lines=len(absence), binary_sha256=row['binary_sha256']))

        if len(result['runs']) != 15 or len(result['gates']) != 12:
            raise ValueError('Run or numerical-gate count differs')
        for name, gate in result['gates'].items():
            report = gate['report']
            if not gate['passed'] or gate['failures'] or gate['exact_observed']:
                raise ValueError('Numerical gate differs: '+name)
            if any(report[key] for key in ('header_differences','topology_differences',
                                           'threshold_max_abs','threshold_nonfinite')):
                raise ValueError('Checkpoint header/topology/threshold differs: '+name)

        logging = result['decisions']['logging']
        loud = [row['median_seconds'] for row in runs[3:11] if row['mode']=='loud']
        quiet = [row['median_seconds'] for row in runs[3:11] if row['mode']=='quiet']
        loud_mean, quiet_mean = statistics.mean(loud), statistics.mean(quiet)
        equal(loud_mean, logging['loud_mean_seconds'])
        equal(quiet_mean, logging['quiet_mean_seconds'])
        equal(loud_mean-quiet_mean, logging['saving_seconds'])
        pairs = [runs[i]['median_seconds']-runs[i+1]['median_seconds'] for i in (3,7)]
        pairs += [runs[i+3]['median_seconds']-runs[i+2]['median_seconds'] for i in (3,7)]
        paired = [pairs[0],pairs[2],pairs[1],pairs[3]]
        for actual, expected in zip(paired, logging['paired_savings_seconds']):
            equal(actual, expected)
        if logging['independent_runs_per_mode'] != 4:
            raise ValueError('Logging run count differs')

        residual = result['decisions']['residual']
        main_mean = statistics.mean([runs[11]['median_seconds'],runs[14]['median_seconds']])
        candidate_mean = statistics.mean([runs[12]['median_seconds'],runs[13]['median_seconds']])
        equal(main_mean, residual['main_mean_seconds'])
        equal(candidate_mean, residual['quiet_candidate_mean_seconds'])
        equal(candidate_mean-main_mean, residual['candidate_minus_main_seconds'])
        equal(candidate_mean/main_mean, residual['candidate_to_main_ratio'])
        for actual, expected in zip([runs[12]['median_seconds']/runs[11]['median_seconds']-1,
                                     runs[13]['median_seconds']/runs[14]['median_seconds']-1],
                                    residual['paired_regression_fractions']):
            equal(actual, expected)

    with path.open('rb') as source:
        digest = hashlib.file_digest(source,'sha256').hexdigest()
    return dict(status='verified', job_id=result['slurm_job_id'],review_archive=str(path),
                review_archive_sha256=digest,verified_files=len(manifest),host=result['host'],
                runs=runs,numerical_gates=len(result['gates']),
                logging=dict(loud_mean_seconds=loud_mean,quiet_mean_seconds=quiet_mean,
                             loud_minus_quiet_seconds=loud_mean-quiet_mean,
                             paired_savings_seconds=paired,independent_runs_per_mode=4),
                residual=dict(main_mean_seconds=main_mean,quiet_candidate_mean_seconds=candidate_mean,
                              candidate_minus_main_seconds=candidate_mean-main_mean,
                              candidate_to_main_ratio=candidate_mean/main_mean,
                              independent_runs_per_version=2))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_250.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]),indent=2))
