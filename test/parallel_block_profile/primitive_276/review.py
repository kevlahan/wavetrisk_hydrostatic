"""Independently verify Stage276 primitive evidence and control timing results."""
import json
import math
from pathlib import Path
import re
import sys
import tarfile
import tempfile
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))
from analyze import analyze
from common import CHECKPOINTS, J5_HASHES, MAIN, BASELINE, ORDER, VALIDATIONS, sha, summarize
from control import verify_archive
from placement_200 import verify
from protocol_229 import ordinary
from validation_230 import storage_evidence
from evidence import publication_evidence
from checkpoint_gate import accept
from contract import evidence as contract_evidence


def same_derived_report(left, right):
    """Permit only float rounding in derived reports, never raw evidence.

    Python 3.12 changed built-in float summation. Cluster 3.11 and local
    Python can therefore differ by a few ulps even with identical raw rows.
    Hashes, keys, ordering, integer counts and all other values remain exact.
    """
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(same_derived_report(left[k], right[k]) for k in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(same_derived_report(a, b) for a, b in zip(left, right))
    if isinstance(left, float) and isinstance(right, float):
        return math.isfinite(left) and math.isfinite(right) and abs(left-right) <= 8*max(math.ulp(left),math.ulp(right))
    return type(left) is type(right) and left == right


def read_runtime(root, result):
    from accounting import parse_accounting
    account=json.loads((root/'accounting.json').read_text())
    parsed=parse_accounting(account['stdout'],str(result['job_id'])) if account['returncode']==0 else None
    if parsed != account['job']:raise ValueError('Accounting report differs from raw output')
    if result.get('elapsed_seconds',0)<=0 or not result.get('started_utc') or not result.get('ended_utc'):
        raise ValueError('Missing campaign elapsed time')
    return dict(slurm=parsed,campaign_started_utc=result['started_utc'],campaign_ended_utc=result['ended_utc'],
                campaign_elapsed_seconds=result['elapsed_seconds'],
                limitation='Slurm elapsed includes batch setup; runner elapsed begins before preflight. Queue and collection/transfer excluded. Accounting can be unavailable; no timestamp inference.')


def review(archive):
    archive = Path(archive).resolve(strict=True)
    count = verify_archive(archive)
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(archive) as stream:
            stream.extractall(tmp, filter='data')
        root = Path(tmp) / 'review'
        out = root / 'results'
        result = json.loads((out / 'results.json').read_text())
        if result['status'] != 'complete':
            raise ValueError('Diagnostic archive: ' + result['status'] + ': ' + str(result.get('error')))
        identity = json.loads((root / 'kit/IDENTITY.json').read_text())
        if (identity['stage'] != 276 or identity.get('revision') != 'primitive-prototype' or result['identity'] != identity or identity['main_commit'] != MAIN or
                identity['stage266_commit'] != BASELINE):
            raise ValueError('Source identity differs')
        if (sha(root / 'kit/files-sha256.json') != result['kit_manifest_sha256'] or
                result['input_hashes'] != J5_HASHES or sha(out / 'input/simple.in') != J5_HASHES['simple.in']):
            raise ValueError('Kit or input seal differs')
        notes = json.loads((out / 'builds.json').read_text())
        if (notes != result['builds'] or notes['status'] != 'passed' or notes['level'] != 'J5' or
                set(notes['builds']) != {'main','stage266','candidate','checked'} or len(notes['physics_objects']) != 23):
            raise ValueError('Build evidence incomplete')
        source_manifest = root / 'kit/source-files.json'
        if sha(source_manifest) != notes['source_manifest_sha256']:
            raise ValueError('Source manifest differs')
        from common import source_scope
        source_scope(root / 'kit')
        for kind, record in notes['builds'].items():
            if (record['physics_objects'] != notes['physics_objects'] or
                    sha(out / ('build-' + kind + '.log')) != record['log_sha256'] or
                    '-j1' not in record['command'] or ('DEBUG=check' in record['command']) != (kind == 'checked')):
                raise ValueError('Build command or shared physics differs')
        if result.get('repeat_used') or len(result['schedule']) != len(VALIDATIONS)+len(ORDER) or set(result['schedule']) != set(result['runs']):
            raise ValueError('Require the complete validation and control schedule')
        validations = [(result['runs'][name]['ranks'],result['runs'][name]['kind'],result['runs'][name]['oracle'])
                       for name in result['schedule'] if result['runs'][name]['phase']=='validation']
        if validations != list(VALIDATIONS):raise ValueError('Required validation schedule differs')
        timing = []
        anchors = {}
        states = None
        for name in result['schedule']:
            row = result['runs'][name]
            directory = out / name
            log = (directory / 'run.log').read_text()
            counts = re.findall(r'number of processors\s*=\s*(\d+)', log)
            if not counts or any(int(n) != row['ranks'] for n in counts):
                raise ValueError('Solver rank count differs: ' + name)
            if sha(directory / 'simple.in') != row['input_sha256'] or row['input_sha256'] != J5_HASHES['simple.in']:
                raise ValueError('Run input differs')
            analysis = analyze(directory / 'run.log')
            raw_states = [step['state'] for step in analysis['step_records']]
            if states is None:
                states = raw_states
            if (not analysis['completed'] or len(raw_states) != 12 or raw_states != states or
                    row['states'] != raw_states or log.count('Remapping vertical coordinates') != 2 or
                    re.findall(r'Restarting from checkpoint\s+(\d+)', log) != ['3', '4']):
                raise ValueError('Science/reload/remap evidence differs: ' + name)
            case = json.loads((out / f'placement-{row["ranks"]}.json').read_text())
            records = [json.loads((directory / f'rank-{rank}.json').read_text()) for rank in range(row['ranks'])]
            verify(records, case)
            if (row['returncode'] != 0 or row['binary_sha256'] != notes['builds'][row['kind']]['binary_sha256'] or
                    row['checkpoint_sha256'] != CHECKPOINTS[row['family']] or sha(directory / 'run.log') != row['log_sha256']):
                raise ValueError('Run identity differs: ' + name)
            numerical = json.loads((directory / 'numerical.json').read_text())
            if (numerical != result['gates'][name] or not numerical['passed'] or
                    accept(numerical['report'], exact=numerical['exact_required'])):
                raise ValueError('Numerical comparison failed')
            if row['family'] == 'stage266':
                exact = json.loads((directory / 'block-exact.json').read_text())
                if not exact['passed'] or accept(exact['report'], exact=True):
                    raise ValueError('Exact block comparison failed')
            vtk = json.loads((directory / 'vtk-members.json').read_text())
            if len(vtk) != 62 or sha(directory / 'vtk-members.json') != row['vtk_manifest_sha256']:
                raise ValueError('VTK evidence incomplete')
            key = (row['ranks'], row['family'])
            if key in anchors and vtk != anchors[key]:
                raise ValueError('Same-count VTK hashes differ')
            anchors[key] = vtk
            expected_catalogue=row['phase']=='validation' and row['kind'] in ('candidate','checked')
            if row.get('catalogue')!=expected_catalogue:raise ValueError('Contract diagnostic identity differs')
            if not same_derived_report(contract_evidence(directory,row['ranks'],expected_catalogue),row['contract']):
                raise ValueError('Raw ownership contract evidence differs')
            if publication_evidence(log, 'main' if row['family']=='main' else 'candidate', row['oracle']) != row['publication']:
                raise ValueError('Direct publication evidence differs')
            storage = storage_evidence(log, row['ranks'], 'baseline' if row['family'] == 'main' else 'candidate', row['oracle'])
            # JSON serializes census tuples as lists; compare the persisted form.
            if json.loads(json.dumps(storage)) != row['storage']:
                raise ValueError('Storage policy differs')
            if row['family'] == 'stage266':
                for marker in ('Stage203 inverse workspace = on', 'Stage203 scalar fields = on'):
                    if marker not in log:
                        raise ValueError('Missing resident inverse selector marker')
                if row['ranks'] > 160 and f'Resident block execution: MPI tasks = {row["ranks"]}, geometric Domains = 160' not in log:
                    raise ValueError('Missing zero-Domain-rank marker')
            if row['oracle']:
                for marker in ('dynamics', 'adaptation', 'vertical-remap'):
                    if f'Block {marker} validation oracle enabled' not in log:
                        raise ValueError('Missing oracle evidence')
            if 'Parallel-block profiling summary (inclusive):' in log:
                raise ValueError('Unexpected profiling')
            if list(directory.glob('cost-rank-*.txt')):raise ValueError('Audit output in uninstrumented run')
            values = ordinary(log, analysis, (directory / 'test_log').read_text())
            if values != row['ordinary'] or [x['index'] for x in values['excluded']] != [0, 7, 8]:
                raise ValueError('Ordinary timing differs')
            if row['phase'] == 'timing':
                if row['oracle'] or row['kind'] == 'checked' or (row['ranks'],row['kind']) not in ORDER:
                    raise ValueError('Invalid timing sample')
                timing.append(dict(row, ordinary=values))
        resident_gate=all(row['contract']['resident_gate_passed'] for row in result['runs'].values() if row['catalogue'])
        if result['prototype_residency_passed']!=resident_gate:raise ValueError('Prototype residency summary differs')
        comparison = summarize(timing)
        if not same_derived_report(comparison, result['comparison']):
            raise ValueError('Timing comparison differs')
        return dict(status='verified', archive_sha256=sha(archive), verified_files=count,
                    identity=identity, prototype_residency_passed=resident_gate, comparison=comparison, runtime=read_runtime(root,result),
                    contracts={name:row['contract'] for name,row in result['runs'].items() if row['kind'] in ('candidate','checked')},
                    derived_report_float_tolerance='Eight ulps only; raw member hashes and all scientific gates remain exact.',
                    limitation='Binary/checkpoint/VTK payloads omitted; cluster-computed hashes and numerical reports verified, not rerun.')


if __name__ == '__main__':
    print(json.dumps(review(sys.argv[1]), indent=2))
