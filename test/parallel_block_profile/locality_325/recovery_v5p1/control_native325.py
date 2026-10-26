"""One new verification allocation and bounded collection; no automatic retry."""
from pathlib import Path
import getpass
import hashlib
import io
import json
import os
import selectors
import signal
import stat
import subprocess
import sys
import tarfile
import tempfile
import time

MAX_TOTAL = 32 * 1024**2
MAX_MEMBER = 2 * 1024**2
MAX_TEST_LOG = 2 * 1024**2
PREFIX = 'rk325-native-v5p1'
TERMINAL = {'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'DEADLINE', 'FAILED', 'NODE_FAIL',
            'OUT_OF_MEMORY', 'PREEMPTED', 'REVOKED', 'SPECIAL_EXIT', 'TIMEOUT'}
ACTIVE = {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED', 'RESIZING',
          'REQUEUED', 'REQUEUE_FED', 'REQUEUE_HOLD', 'RESV_DEL_HOLD', 'SIGNALING',
          'STAGE_OUT', 'STOPPED'}


def require(ok, detail):
    if not ok:
        raise ValueError('native.control325.' + detail)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def safe(root, name):
    rel = Path(name)
    require(type(name) is str and not rel.is_absolute() and '..' not in rel.parts and
            str(rel) == name, 'unsafe_relative_path')
    root = Path(root).resolve(strict=True)
    path = root / rel
    require(path.parent.resolve(strict=True).is_relative_to(root) and not path.is_symlink(),
            'unsafe_input_path')
    return path


def package_check(root):
    root = Path(root).resolve(strict=True)
    manifest_path = root / 'native-files-sha256.json'
    require(manifest_path.is_file() and not manifest_path.is_symlink(), 'package_manifest_required')
    manifest = json.loads(manifest_path.read_text())
    require(type(manifest) is dict and bool(manifest), 'package_manifest_schema')
    actual = set()
    for base, directories, files in os.walk(root, followlinks=False):
        for name in directories + files:
            path = Path(base) / name
            info = path.lstat()
            require(not stat.S_ISLNK(info.st_mode) and
                    (stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode)), 'package_entry_type')
            if stat.S_ISREG(info.st_mode) and path != manifest_path:
                actual.add(str(path.relative_to(root)))
    require(actual == set(manifest), 'package_inventory')
    for name, digest in manifest.items():
        path = safe(root, name)
        require(path.is_file() and type(digest) is str and len(digest) == 64 and
                all(c in '0123456789abcdef' for c in digest) and sha(path) == digest,
                'package_bytes:' + name)
    return sha(manifest_path)


def save(path, value):
    path = Path(path)
    data = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    require(len(data) <= MAX_MEMBER and not path.is_symlink(), 'state_bound_or_link')
    temp = path.with_name(path.name + '.tmp')
    require(not temp.exists() and not temp.is_symlink(), 'unreviewed_state_temporary')
    descriptor = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'wb') as output:
        output.write(data); output.flush(); os.fsync(output.fileno())
    os.replace(temp, path)
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def queued(job):
    require(not isinstance(job, bool) and str(job).isdigit() and int(job) > 0, 'scheduler_job_id')
    job = str(job)
    q = subprocess.run(['squeue', '-h', '-u', getpass.getuser(), '-o', '%A'],
                       capture_output=True, text=True)
    if q.returncode == 0 and job in q.stdout.split():
        return True
    account = subprocess.run(['sacct', '-n', '-P', '-X', '-j', job, '-o', 'JobIDRaw,State'],
                             capture_output=True, text=True)
    states = []
    for line in account.stdout.splitlines():
        fields = line.strip().split('|')
        if len(fields) > 1 and fields[0] == job and fields[1].split():
            states.append(fields[1].split()[0].rstrip('+'))
    if account.returncode == 0 and states and all(value in TERMINAL for value in states):
        return False
    if account.returncode == 0 and any(value in ACTIVE for value in states):
        return True
    raise RuntimeError('Scheduler state unknown; stop without collection/retry. squeue=' +
                       q.stderr + '; sacct=' + account.stderr + ' ' + account.stdout)


def state_path():
    return Path.home() / (PREFIX + '-latest.json')


def record(state):
    save(state_path(), state)
    if state.get('work'):
        save(Path(state['work']) / 'submission.json', state)


def run_preflight(kit, log, *, timeout=300):
    """Bound log bytes before writes and terminate only our own test group."""
    kit, log = Path(kit), Path(log)
    require(not log.exists() and not log.is_symlink(), 'fresh_preflight_log')
    command = [sys.executable, '-B', '-m', 'unittest', 'discover', '-v']
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', BLIS_NUM_THREADS='1')
    process = subprocess.Popen(command, cwd=kit, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, start_new_session=True, env=environment)
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    started, count = time.monotonic(), 0
    try:
        with log.open('xb') as output:
            while selector.get_map():
                require(time.monotonic() - started < timeout, 'preflight_timeout')
                for key, _ in selector.select(timeout=.25):
                    data = os.read(key.fileobj.fileno(), 16384)
                    if not data:
                        selector.unregister(key.fileobj)
                        continue
                    require(count + len(data) <= MAX_TEST_LOG, 'preflight_log_cap')
                    output.write(data)
                    count += len(data)
            process.wait(timeout=5)
            output.flush(); os.fsync(output.fileno())
        require(process.returncode == 0, 'preflight_tests_failed')
        return dict(passed=True, command=command, returncode=process.returncode,
                    log_bytes=count, log_sha256=sha(log), elapsed_seconds=time.monotonic()-started)
    finally:
        selector.close()
        active_error = sys.exc_info()[1]
        cleanup_errors = []
        def send(signum):
            try:
                os.killpg(process.pid, signum)
            except ProcessLookupError:
                pass
            except OSError as error:
                cleanup_errors.append(repr(error))
        send(signal.SIGTERM)
        # A test parent can already have exited while an ignoring grandchild
        # keeps its stdout pipe and process group alive. Escalate the group
        # independently of parent.poll(), then reap the parent without masking
        # the original test/cap/timeout error if cleanup itself has trouble.
        try:
            process.wait(timeout=.2)
        except subprocess.TimeoutExpired:
            pass
        time.sleep(.05)
        send(signal.SIGKILL)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired as cleanup_error:
            if active_error is None:
                raise RuntimeError('Preflight process group did not stop; no submission') from cleanup_error
            if hasattr(active_error, 'add_note'):
                active_error.add_note('Preflight cleanup did not reap its parent within five seconds')
        finally:
            process.stdout.close()
        if cleanup_errors:
            detail = 'Preflight group cleanup errors: ' + '; '.join(cleanup_errors)
            if active_error is None:
                raise RuntimeError(detail)
            if hasattr(active_error, 'add_note'):
                active_error.add_note(detail)


def submit(kit):
    kit = Path(kit).resolve(strict=True)
    package_sha = package_check(kit)
    for job in (55569, 55571, 55572, 55573, 55574, 55575):
        require(not queued(job), 'predecessor_still_active:' + str(job))
    lock = Path.home() / ('.' + PREFIX + '-launch.lock')
    lock.mkdir(mode=0o700)
    try:
        require(not state_path().exists(), 'state_exists_no_automatic_retry')
        intent = Path.home() / (PREFIX + '-launch-intent.json')
        require(not intent.exists() and not intent.is_symlink(), 'submission_intent_exists')
        work = Path(tempfile.mkdtemp(prefix=PREFIX + '-', dir=Path.home()))
        state = dict(schema=1, stage=325, revision='native-qualification-v5p1', source_job=55569,
                     predecessor_job=55575, status='preflight', job=None, kit=str(kit),
                     package_manifest_sha256=package_sha, work=str(work), out=str(work / 'out'),
                     log=str(work / 'slurm.log'), preflight_log=str(work / 'preflight-tests.log'),
                     walltime_seconds=1800, builds=0, solver_runs=0, performance_runs=0,
                     predecessor_elapsed_seconds=0, original_budget_used_seconds=2656,
                     new_budget_remaining_seconds=1800,
                     automatic_retry=False, independent_review_required=True)
        save(intent, dict(status='preflight_before_single_submission', state=state))
        record(state)
        try:
            state['preflight'] = run_preflight(kit, Path(state['preflight_log']))
            require(package_check(kit) == package_sha, 'package_changed_during_tests')
        except BaseException as error:
            state.update(status='preflight_failed', error=repr(error))
            record(state)
            raise
        state['status'] = 'submitting'
        record(state)
        command = ['sbatch', '--parsable', '--chdir=' + str(work), '--output=' + state['log'],
                   str(kit / 'verification.sbatch'), str(kit), state['out']]
        save(work / 'sbatch-intent.json', dict(command=command, status='submitting'))
        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except BaseException as error:
            state.update(status='submission_unknown', error=repr(error))
            record(state)
            raise
        save(work / 'sbatch-response.json', dict(command=command, returncode=result.returncode,
                                                stdout=result.stdout, stderr=result.stderr))
        answer = result.stdout.strip().split(';')[0]
        if result.returncode != 0 or not answer.isdigit() or int(answer) <= 0:
            state.update(status='submission_unknown', error='No confirmed job ID; inspect scheduler, never retry automatically')
            record(state)
            raise RuntimeError(state['error'])
        state.update(status='submitted', job=int(answer))
        record(state)
        print(json.dumps(state, indent=2), flush=True)
        return state
    finally:
        lock.rmdir()


def status(state):
    result = dict(status=state['status'], job=state.get('job'))
    if state.get('job'):
        result['active'] = queued(state['job'])
    for path in (Path(state['out']) / 'supervisor-report.json',
                 Path(state['out']) / 'state.attempt.json'):
        if path.is_file() and not path.is_symlink():
            require(path.stat().st_size <= MAX_MEMBER, 'status_report_bound')
            obj = json.loads(path.read_text())
            result[path.name] = {key: obj[key] for key in
                                 ('status', 'phase', 'elapsed_seconds', 'signal', 'error',
                                  'native_gate_passed', 'native_qualification') if key in obj}
    print(json.dumps(result, indent=2), flush=True)
    return result


def selected_files(state):
    kit, work = Path(state['kit']), Path(state['work'])
    package_check(kit)
    manifest = json.loads((kit / 'native-files-sha256.json').read_text())
    files, excluded, total = [], [], 0
    def add(path, name):
        nonlocal total
        stat = path.stat()
        require(not path.is_symlink() and path.is_file() and stat.st_size <= MAX_MEMBER,
                'collection_member_bound:' + name)
        total += stat.st_size
        require(total <= MAX_TOTAL - 2 * MAX_MEMBER, 'collection_total_bound')
        files.append((path, name, stat.st_size, stat.st_mtime_ns))
    for name in sorted(set(manifest) | {'native-files-sha256.json'}):
        add(safe(kit, name), 'kit/' + name)
    for path in sorted(work.rglob('*')):
        relative = path.relative_to(work)
        if not path.is_file() and not path.is_symlink():
            continue
        omitted = ('__pycache__' in relative.parts or 'cache-build' in relative.parts or
                   'work-input' in relative.parts or path.suffix not in ('.json', '.log') or
                   path.name == 'work325-rank-vectors.json' or
                   path.name.startswith(('execution316-', 'test_checkpoint_')))
        if path.is_symlink() or omitted:
            excluded.append(str(relative))
            continue
        add(path, 'work/' + str(relative))
    return files, excluded, total


def verify_archive(path):
    with tarfile.open(path) as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        require(len(names) == len(set(names)) and all(member.isfile() and
                not Path(member.name).is_absolute() and '..' not in Path(member.name).parts
                for member in members), 'review_member_safety')
        manifest = json.load(archive.extractfile('review/files-sha256.json'))
        require(set(names) == {'review/files-sha256.json'} | {'review/' + n for n in manifest},
                'review_inventory')
        require(sum(member.size for member in members) <= MAX_TOTAL, 'review_logical_cap')
        for name, digest in manifest.items():
            require(hashlib.file_digest(archive.extractfile('review/' + name), 'sha256').hexdigest() == digest,
                    'review_file_hash:' + name)
    return len(manifest)


def collect(state):
    if state.get('job'):
        require(not queued(state['job']), 'terminal_job_required')
        accounting = subprocess.run(['sacct', '-n', '-P', '-X', '-j', str(state['job']), '-o',
                                     'JobIDRaw,State,ExitCode,Start,End,Elapsed,ElapsedRaw,NodeList'],
                                    capture_output=True, text=True)
        require(accounting.returncode == 0, 'accounting_unavailable')
        rows = [line.split('|') for line in accounting.stdout.splitlines() if line.strip()]
        require(any(row[0] == str(state['job']) and row[1].split()[0].rstrip('+') in TERMINAL
                    for row in rows if len(row) >= 8), 'accounting_terminal_receipt_required')
        account = dict(command=accounting.args, returncode=accounting.returncode,
                       stdout=accounting.stdout, stderr=accounting.stderr)
        label = str(state['job'])
    else:
        require(state['status'] == 'preflight_failed', 'confirmed_submission_or_preflight_failure_required')
        account, label = dict(status='not_submitted'), 'preflight'
    destination = Path.home() / ('rk325-native-v5p1-review-' + label + '.tar.gz')
    if destination.exists():
        with tarfile.open(destination) as archive:
            prior = json.load(archive.extractfile('review/submission.json'))
            require(prior.get('job') == state.get('job') and
                    prior.get('package_manifest_sha256') == state.get('package_manifest_sha256'),
                    'existing_review_identity')
        print('Verified', verify_archive(destination), 'files', destination, 'SHA256', sha(destination), flush=True)
        return destination
    files, excluded, size = selected_files(state)
    generated = {'submission.json': state, 'accounting.json': account,
                 'collection.json': dict(selected_logical_bytes=size, maximum_logical_bytes=MAX_TOTAL,
                                         excluded_paths=excluded, excluded_contents_read=False,
                                         policy='Sealed code and bounded JSON/logs; SQLite, captures, checkpoints, TSV and vectors remain native.')}
    temp = Path(str(destination) + '.tmp')
    require(not temp.exists() and not temp.is_symlink(), 'collection_temporary_requires_review')
    manifest, written = {}, 0
    with tarfile.open(temp, 'x:gz') as archive:
        def add(name, data):
            nonlocal written
            require(len(data) <= MAX_MEMBER, 'generated_or_selected_member_bound')
            require(written + len(data) <= MAX_TOTAL - MAX_MEMBER, 'archive_write_cap')
            info = tarfile.TarInfo('review/' + name)
            info.size, info.mode = len(data), 0o600
            archive.addfile(info, io.BytesIO(data))
            manifest[name] = hashlib.sha256(data).hexdigest()
            written += len(data)
        for name, value in generated.items():
            add(name, (json.dumps(value, sort_keys=True, indent=2) + '\n').encode())
        for path, name, expected_size, expected_mtime in files:
            data = path.read_bytes()
            stat = path.stat()
            require(not path.is_symlink() and stat.st_size == len(data) == expected_size and
                    stat.st_mtime_ns == expected_mtime, 'collection_input_changed:' + name)
            add(name, data)
        data = (json.dumps(manifest, sort_keys=True, indent=2) + '\n').encode()
        require(len(data) <= MAX_MEMBER and written + len(data) <= MAX_TOTAL,
                'archive_manifest_cap')
        info = tarfile.TarInfo('review/files-sha256.json')
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    verify_archive(temp)
    with temp.open('rb') as raw:
        os.fsync(raw.fileno())
    os.rename(temp, destination)
    print('Verified', verify_archive(destination), 'files', destination, 'SHA256', sha(destination), flush=True)
    return destination


if __name__ == '__main__':
    require(len(sys.argv) == 3 and sys.argv[1] in ('submit', 'status', 'collect'),
            'usage: control_native325.py submit|status|collect KIT')
    if sys.argv[1] == 'submit':
        submit(sys.argv[2])
    else:
        state = json.loads(state_path().read_text())
        require(Path(state['kit']).resolve() == Path(sys.argv[2]).resolve(), 'controller_kit_state_mismatch')
        (status if sys.argv[1] == 'status' else collect)(state)
