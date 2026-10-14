"""One-command launch/status/collection, including failed-job diagnostics."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
from common import preflight, save, sha, verify
from accounting import capture


def state_path():
    return Path.home() / 'rk278-latest.json'


def queued(job):
    query = subprocess.run(['squeue', '-h', '-u', os.environ['USER'], '-o', '%A'],
                           capture_output=True, text=True, check=True)
    return str(job) in query.stdout.split()


def submit(kit, repeat=False):
    kit = Path(kit).resolve(strict=True)
    if repeat: raise ValueError('Repeats require review of the initial campaign and a separate bounded schedule')
    home = Path.home()
    fixture = home / 'stage187-fixture-20260911-145307'
    preflight(kit, fixture)
    if not shutil.which('sbatch'):
        raise ValueError('Run launch.sh on bbserv with Slurm available')
    lock = home / '.rk278-launch.lock'
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError('Another Stage278 launch is in progress: ' + str(lock)) from None
    try:
        latest = state_path()
        if latest.exists():
            prior = json.loads(latest.read_text())
            archive = home / f'rk278-review-{prior["job"]}.tar.gz'
            if not archive.is_file():
                raise ValueError(f'Job {prior["job"]} already submitted; run status.sh and collect.sh before another launch')
            verify_archive(archive)
        work=Path(tempfile.mkdtemp(prefix='rk278-attribution-',dir=home))
        out=work/'results'
        response = subprocess.check_output(['sbatch', '--parsable', '--chdir=' + str(home),
            '--output=' + str(home / 'rk278-%j.log'), str(kit / 'timing.sbatch'),
            str(kit), str(fixture), str(out)], text=True).strip().split(';')[0]
        if not response.isdigit():
            raise ValueError('Unexpected sbatch response: ' + response)
        state = dict(phase='initial', job=int(response), kit=str(kit), out=str(out), fixture=str(fixture),
                     log=str(home / f'rk278-{response}.log'))
        save(home / f'rk278-job-{response}.json', state)
        save(latest, state)
        print(f'Submitted attribution-only campaign on nk/bb13: {response}')
        print('Log: ' + state['log'])
        print('After completion: bash ' + str(kit / 'collect.sh'))
    finally:
        lock.rmdir()


def status(state):
    job = state['job']
    print(f'Job {job}: ' + ('queued/running' if queued(job) else 'not in active queue'))
    accounting = subprocess.run(['sacct', '-j', str(job), '-X',
                                '-o', 'JobID,JobName%28,State,ExitCode,Elapsed,NodeList'],
                               capture_output=True, text=True)
    if accounting.returncode == 0:
        print(accounting.stdout.rstrip())
    else:
        print('Accounting unavailable; using saved log and results.')
    out = Path(state['out'])
    if (out / 'results.json').exists():
        result = json.loads((out / 'results.json').read_text())
        print('status:', result['status'], 'active:', result.get('active'), 'error:', result.get('error'))
    log = Path(state['log'])
    print('Log: ' + str(log))
    if log.exists():
        print('\n'.join(log.read_text(errors='replace').splitlines()[-35:]))


def verify_archive(archive):
    with tarfile.open(archive) as stream:
        manifest = json.load(stream.extractfile('review/files-sha256.json'))
        for name, expected in manifest.items():
            member = stream.extractfile('review/' + name)
            if member is None or __import__('hashlib').file_digest(member, 'sha256').hexdigest() != expected:
                raise ValueError('Review member differs: ' + name)
    return len(manifest)


def collect(state):
    if queued(state['job']):
        raise ValueError('Job is still active; collect after it ends')
    home = Path.home()
    archive = home / f'rk278-review-{state["job"]}.tar.gz'
    if not archive.exists():
        with tempfile.TemporaryDirectory(prefix='.rk278-collect-', dir=home) as tmp:
            review = Path(tmp) / 'review'
            review.mkdir()
            save(review / 'submission.json', state)
            save(review / 'accounting.json', capture(state['job']))
            log = Path(state['log'])
            if log.exists():
                shutil.copy2(log, review / 'slurm.log')
            kit = Path(state['kit'])
            (review / 'kit').mkdir()
            for name in ('IDENTITY.json', 'source-files.json', 'files-sha256.json', 'LOCAL_VALIDATION.json', 'TEST_VALIDATION.json', 'ENTRY_GATE.json', 'overlay-files.json', 'campaign-spec.json', 'parent-source-files.json', 'diagnostic.py', 'attribution.py', 'review.py', 'run.py', 'build.py', 'common.py'):
                if (kit / name).exists():
                    shutil.copy2(kit / name, review / 'kit' / name)
            out = Path(state['out'])
            if out.exists():
                for path in sorted(out.rglob('*')):
                    if (not path.is_file() or path.is_symlink() or 'sources' in path.relative_to(out).parts):
                        continue
                    if path.suffix not in ('.json', '.log', '.txt', '.resources') and path.name not in ('test_log', 'simple.in'):
                        continue
                    destination = review / 'results' / path.relative_to(out)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(path, destination)
            files = {str(p.relative_to(review)): sha(p) for p in sorted(review.rglob('*')) if p.is_file()}
            save(review / 'files-sha256.json', files)
            temporary = archive.with_suffix('.tmp')
            with tarfile.open(temporary, 'w:gz', format=tarfile.PAX_FORMAT) as stream:
                stream.add(review, arcname='review')
            verify_archive(temporary)
            temporary.rename(archive)
    print(f'PASS: {verify_archive(archive)} review files; SHA256 {sha(archive)}')
    subprocess.run(['scp', str(archive), 'kevlahan@murier.math.mcmaster.ca:~/Desktop/'], check=True)
    print('Transferred review to murier Desktop: ' + archive.name)


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] not in ('submit', 'status', 'collect'):
        raise SystemExit('Usage: control.py submit|status|collect KIT')
    action, kit = sys.argv[1:]
    if action=='submit':
        submit(kit)
    else:
        state = json.loads(state_path().read_text())
        if Path(state['kit']).resolve()!=Path(kit).resolve():
            raise ValueError('Latest job belongs to another kit: '+state['kit'])
        (status if action == 'status' else collect)(state)
