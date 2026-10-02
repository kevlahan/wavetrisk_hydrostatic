"""Independently review the Stage254 partial return, including its usable probe log."""
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
            'instrumented': 'f0d5daac2da9dfaa4dc8027e8ab8ad1147c06ca088818aee097f605bf334ca22'}
OBJECT_MAP_SHA = '9717172ac48acbe10744240b539d1f562f65bc4d230a6fc2228f64945857f031'
PHASES = ('complete timestep', 'native velocity source', 'native velocity gradient')
OWNER_PHASES = ('execute', 'resident fetch', 'internal read', 'read MPI',
                'producer loop', 'delivery read', 'delivery MPI')
PROFILE_ROW = re.compile(r'^\s*(.+?):\s*(\d+)\s*\|\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*\|')
OWNER_ROW = re.compile(r'^Owner program (.+?) seconds min/avg/max =\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)')
DEPTH_ROW = re.compile(r'^Owner program depth (\d+) calls =\s*(\d+) read seconds avg =\s*([\d.Ee+-]+)')


def review(path):
    path = Path(path).resolve(strict=True)
    with tarfile.open(path) as tar:
        listed = [member for member in tar if member.isfile()]
        members = {member.name: member for member in listed}
        if len(members) != len(listed):
            raise ValueError('Duplicate archive member')
        manifest = json.load(tar.extractfile(members['results/REVIEW_SHA256.json']))
        if set(members) != {'results/'+name for name in manifest} | {'results/REVIEW_SHA256.json'}:
            raise ValueError('Archive inventory differs')
        for name, expected in manifest.items():
            if hashlib.sha256(tar.extractfile(members['results/'+name]).read()).hexdigest() != expected:
                raise ValueError('Archive member changed: '+name)

        def read_json(name):
            return json.load(tar.extractfile(members['results/'+name]))

        def read_text(name):
            return tar.extractfile(members['results/'+name]).read().decode()

        result = read_json('results.json')
        expected = [f'J5-n83-{mode}-{index}-{kind}' for mode in ('off', 'profile')
                    for index, kind in enumerate(('baseline', 'instrumented', 'instrumented', 'baseline'), 1)][:6]
        if (result['status'] != 'failed' or result['job_id'] != '55260' or
                result['stage251_job_id'] != '55256' or result['host'] != 'bb13.rhpcs.internal' or
                result['error'] != "ValueError('Incomplete owner timing windows')" or
                result['schedule'] != expected or list(result['runs']) != expected or
                result['active'] != expected[-1]):
            raise ValueError('Campaign identity or failure differs')
        build = result['build']
        if (build['baseline_binary_sha256'] != BINARIES['baseline'] or
                build['instrumented_binary_sha256'] != BINARIES['instrumented'] or
                build['command'][:2] != ['make', '-j1'] or
                'BLOCK_VELOCITY_OWNER=on' not in build['command'] or
                len(result['regenerated_physics_files']) != 23):
            raise ValueError('Build provenance differs')
        objects_sha = hashlib.sha256(json.dumps(build['physics_objects'], sort_keys=True,
                                               separators=(',', ':')).encode()).hexdigest()
        if objects_sha != OBJECT_MAP_SHA:
            raise ValueError('Physics object map differs')
        mapping = {item['rank']: item['cpu'] for item in read_json('placement-83.json')['mapping']}
        if set(mapping) != set(range(83)) or len(set(mapping.values())) != 83:
            raise ValueError('Core placement differs')

        medians = {}
        profile_totals = {}
        owner = None
        for name in expected:
            row = result['runs'][name]
            log = read_text(name+'/run.log')
            if (row['returncode'] or row['ranks'] != 83 or
                    row['binary_sha256'] != BINARIES[row['kind']] or
                    hashlib.sha256(log.encode()).hexdigest() != row['log_sha256']):
                raise ValueError('Binary, exit or log differs: '+name)
            lines = [line for line in log.splitlines() if STEP.match(line)]
            states = [re.sub(r'\s+cpu =.*', '', line.strip()) for line in lines]
            if (len(states) != 12 or
                    hashlib.sha256('\n'.join(states).encode()).hexdigest() != STATE or
                    log.count('Remapping vertical coordinates') != 2):
                raise ValueError('Raw J5 trajectory differs: '+name)
            for rank in range(83):
                record = read_json(f'{name}/rank-{rank}.json')
                if (record['rank'] != rank or record['returncode'] or
                        record['hostname'] != result['host'] or record['affinity'] != [mapping[rank]]):
                    raise ValueError('Rank placement or exit differs: '+name)
            seconds = [float(STEP.match(line)[1]) for line in lines]
            medians[name] = statistics.median(seconds[i] for i in (1,2,3,4,5,6,9,10,11))
            if 'ordinary' in row:
                # The printed CPU field has only two decimal places; the
                # runner's nanosecond timing came from test_log (not archived).
                if (abs(medians[name]-row['ordinary']['median']) > 0.0051 or
                        row['state_sha256'] != STATE or row['checkpoint_sha256'] != CHECKPOINT or
                        row['steps'] != 12 or row['remaps'] != 2 or row['vtk_members'] != 62):
                    raise ValueError('Recorded run gates differ: '+name)
            if row['mode'] == 'off':
                if 'Parallel-block profiling summary (inclusive):' in log or 'Owner program calls execute' in log:
                    raise ValueError('Unprofiled run emitted a timing report: '+name)
                continue
            if log.count('Parallel-block profiling summary (inclusive):') != 2:
                raise ValueError('Profile window count differs: '+name)
            found = {phase: [] for phase in PHASES}
            for line in log.splitlines():
                match = PROFILE_ROW.match(line)
                if match and match[1] in found:
                    found[match[1]].append((int(match[2]), float(match[4])))
            if ([x[0] for x in found['complete timestep']] != [8,4] or
                    [x[0] for x in found['native velocity source']] != [2880,1320] or
                    [x[0] for x in found['native velocity gradient']] != [960,480]):
                raise ValueError('Profile call counts differ: '+name)
            profile_totals[name] = {phase: sum(x[1] for x in values) for phase, values in found.items()}
            if row['kind'] == 'baseline':
                if 'Owner program calls execute' in log:
                    raise ValueError('Owner timer emitted in baseline')
                continue
            calls = [list(map(int, line.split('= ', 1)[1].split())) for line in log.splitlines()
                     if line.startswith('Owner program calls execute/read/delivery = ')]
            if calls != [[2160,10260,2160],[990,4950,990]]:
                raise ValueError('Owner timer call counts differ')
            measured = {phase: [] for phase in OWNER_PHASES}
            depth = {}
            for line in log.splitlines():
                match = OWNER_ROW.match(line)
                if match and match[1] in measured:
                    measured[match[1]].append(tuple(map(float, match.groups()[1:])))
                match = DEPTH_ROW.match(line)
                if match:
                    depth.setdefault(int(match[1]), []).append((int(match[2]),float(match[3])))
            if any(len(values) != 2 for values in measured.values()) or set(depth) != set(range(5)):
                raise ValueError('Owner timing report is incomplete')
            owner = {'calls': calls,
                     'seconds_avg_sum': {phase: sum(v[1] for v in values)
                                         for phase, values in measured.items()},
                     'depth_read_seconds_avg_sum': {str(d): sum(x[1] for x in values)
                                                    for d, values in depth.items()}}

        off = [result['runs'][name]['ordinary']['median'] for name in expected[:4]]
        baseline = statistics.mean((off[0], off[3]))
        instrumented = statistics.mean((off[1], off[2]))
        recorded = result['decisions']['off']
        if (abs(baseline-recorded['baseline_mean_seconds']) > 1e-8 or
                abs(instrumented-recorded['instrumented_mean_seconds']) > 1e-8):
            raise ValueError('Unprofiled A-B-B-A decision differs')
        slurm = read_text('slurm-job.log')
        if 'Incomplete owner timing windows' not in slurm or 'Starting J5-n83-profile-2-instrumented' not in slurm:
            raise ValueError('Slurm failure differs')

    with path.open('rb') as stream:
        archive_sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'status': 'verified_partial', 'job_id': result['job_id'], 'archive_sha256': archive_sha,
            'verified_files': len(manifest), 'completed_runs': 5, 'usable_probe_runs': 1,
            'unprofiled_baseline_seconds': baseline,
            'unprofiled_instrumented_seconds': instrumented,
            'unprofiled_probe_overhead_seconds': instrumented-baseline,
            'profiled_single_run_source_seconds': profile_totals,
            'owner': owner,
            'parser_error': 'Expected all source calls [2880,1320], but measured external-owner calls [2160,990].'}


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: review_cluster_254.py REVIEW_ARCHIVE')
    print(json.dumps(review(sys.argv[1]), indent=2))
