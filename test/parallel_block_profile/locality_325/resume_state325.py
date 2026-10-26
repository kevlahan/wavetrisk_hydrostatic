"""Closed, resumable certificates for *derived* Stage325 verification state.

The runner must freshly hash the protected source, checker, configuration and
capture inventories before opening this store. A progress record never proves a
phase. Only a closed, fsynced directory with an exact certificate can do that.

``reserve`` accounts bytes before external writes; it is not an operating-system
quota. SQLite callers must apply its returned byte allowance as max_page_count,
disable unbudgeted journals, and budget every temporary file/sidecar. The store
checks the whole tree and every actual phase byte again before certification.
"""
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat


PHASES325 = {
    'cache': (),
    'providers-C0': ('cache',),
    'providers-L256': ('cache',),
    'pair': ('cache', 'providers-C0', 'providers-L256'),
    'work': ('pair',),
    'preservation': ('work',),
}
BINDING_KEYS325 = {'source', 'checker', 'config', 'captures'}
_SHA = re.compile(r'[0-9a-f]{64}\Z')
_NAME = re.compile(r'[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}\Z')
_ROOT_FILES = {'identity.json', 'session.json', 'reservations.json',
               'semantic-failure.json', '.lock', '.atomic.tmp'}
_META_SLOTS = 16


class ResumeStateError325(ValueError):
    pass


def require325(ok, detail):
    if not ok:
        raise ResumeStateError325('resume325.' + detail)


def canonical325(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False).encode('ascii') + b'\n'


def sha256325(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _fields(value, fields, label):
    require325(type(value) is dict and set(value) == set(fields), label + '_fields')


def _digest(value):
    return hashlib.sha256(canonical325(value)).hexdigest()


def _safe_absolute(path):
    path = Path(path)
    require325(path.is_absolute() and '..' not in path.parts, 'absolute_clean_path')
    # resolve() alone would silently accept and hide an unsafe parent symlink.
    for component in (path, *path.parents):
        require325(not component.is_symlink(), 'symlink_path')
    return path


def _relative(value):
    require325(type(value) is str and bool(value), 'relative_path')
    path = Path(value)
    require325(not path.is_absolute() and str(path) == value and
               all(_NAME.fullmatch(part) and part not in ('.', '..')
                   for part in path.parts), 'relative_path')
    require325(path.name != 'certificate.json', 'reserved_artifact_path')
    return path


def _fsync_dir(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


class PhaseStore325:
    """Single-writer store; existing state requires ``explicit_resume=True``.

    Binding is exactly four SHA256 strings named source/checker/config/captures.
    ``start`` returns a private staging directory. Reserve a total phase budget
    before writing there, close all artifact writers, then ``commit`` with an
    exact name -> relative path inventory. No writer may remain open at commit.
    """

    def __init__(self, root, binding, *, protected_roots=(), cap_bytes=4 * 1024**3,
                 metadata_limit=1024**2, explicit_resume=False):
        _fields(binding, BINDING_KEYS325, 'binding')
        require325(all(type(x) is str and _SHA.fullmatch(x)
                       for x in binding.values()), 'binding_digest')
        require325(type(cap_bytes) is int and type(metadata_limit) is int and
                   metadata_limit >= 512 and cap_bytes > _META_SLOTS * metadata_limit,
                   'budget_config')
        self.root = _safe_absolute(root)
        self.binding = dict(binding)
        self.binding_sha256 = _digest(binding)
        self.cap_bytes = cap_bytes
        self.metadata_limit = metadata_limit
        self.payload_cap = cap_bytes - _META_SLOTS * metadata_limit
        self.explicit_resume = explicit_resume
        self.protected_roots = tuple(_safe_absolute(p) for p in protected_roots)
        for protected in self.protected_roots:
            require325(not self.root.is_relative_to(protected) and
                       not protected.is_relative_to(self.root), 'protected_root_overlap')
        self._lock = None
        self.certificates = {}
        self._reservations = {}
        self._session = None

    def _writable(self):
        require325(self._lock is not None, 'closed_store')
        require325(not (self.root / 'semantic-failure.json').exists(), 'semantic_failed_latch')

    def available_bytes(self):
        """Unreserved payload space, after all retained payload and metadata headroom."""
        require325(self._lock is not None, 'closed_store')
        return self.payload_cap - sum(self._reservations.values())

    def __enter__(self):
        require325(self._lock is None, 'already_open')
        self.root.mkdir(mode=0o700, parents=False, exist_ok=True)
        root_stat = self.root.stat()
        require325(stat.S_ISDIR(root_stat.st_mode) and root_stat.st_uid == os.getuid()
                   and stat.S_IMODE(root_stat.st_mode) == 0o700, 'private_root')
        lock = self.root / '.lock'
        flags = os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW
        descriptor = os.open(lock, flags, 0o600)
        try:
            require325(stat.S_ISREG(os.fstat(descriptor).st_mode), 'regular_lock')
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ResumeStateError325('resume325.concurrent_writer') from error
            self._lock = descriptor
            self._open(root_stat)
            return self
        except BaseException:
            os.close(descriptor)
            self._lock = None
            raise

    def __exit__(self, kind, value, traceback):
        if self._lock is not None:
            fcntl.flock(self._lock, fcntl.LOCK_UN)
            os.close(self._lock)
            self._lock = None

    def _open(self, root_stat):
        identity = {'schema': 1, 'root': str(self.root), 'device': root_stat.st_dev,
                    'inode': root_stat.st_ino, 'owner': os.getuid(),
                    'binding': self.binding, 'binding_sha256': self.binding_sha256,
                    'cap_bytes': self.cap_bytes, 'metadata_limit': self.metadata_limit}
        existing = (self.root / 'identity.json').exists()
        self._inventory()
        if existing:
            require325(self.explicit_resume is True, 'explicit_resume_required')
            require325(canonical325(self._read('identity.json')) == canonical325(identity),
                       'identity_or_binding_changed')
            # An interrupted atomic metadata update has no certificate authority.
            # Never replace an existing final file with it during resume.
            temp = self.root / '.atomic.tmp'
            if temp.exists():
                require325(temp.is_file() and not temp.is_symlink(), 'unsafe_atomic_temp')
                temp.unlink()
                _fsync_dir(self.root)
            if (self.root / 'semantic-failure.json').exists():
                failure = self._read('semantic-failure.json')
                _fields(failure, ('schema', 'status', 'binding_sha256', 'phase', 'detail'),
                        'semantic_failure')
                require325(type(failure['schema']) is int and failure['schema'] == 1
                           and failure['status'] == 'semantic_failed'
                           and failure['binding_sha256'] == self.binding_sha256,
                           'semantic_failure_schema')
                raise ResumeStateError325('resume325.semantic_failed_latch: ' +
                                          str(failure['detail']))
            reservations = self._read('reservations.json')
            _fields(reservations, ('schema', 'binding_sha256', 'allocations'), 'reservations')
            require325(type(reservations['schema']) is int and reservations['schema'] == 1 and
                       reservations['binding_sha256'] == self.binding_sha256 and
                       type(reservations['allocations']) is dict and
                       set(reservations['allocations']) <= set(PHASES325) and
                       all(type(x) is int and x >= 0
                           for x in reservations['allocations'].values()), 'reservations_schema')
            self._reservations = reservations['allocations']
            self._session = self._read('session.json')
            self._check_session()
        else:
            require325(set(p.name for p in self.root.iterdir()) == {'.lock'},
                       'unidentified_existing_state')
            (self.root / 'phases').mkdir(mode=0o700)
            (self.root / 'partial').mkdir(mode=0o700)
            self._atomic('identity.json', identity)
            self._session = {'schema': 1, 'binding_sha256': self.binding_sha256,
                             'status': 'idle', 'phase': None, 'signal': None}
            self._save_reservations()
            self._save_session()
        self.validate_all()

    def _check_session(self):
        value = self._session
        _fields(value, ('schema', 'binding_sha256', 'status', 'phase', 'signal'), 'session')
        require325(type(value['schema']) is int and value['schema'] == 1
                   and value['binding_sha256'] == self.binding_sha256
                   and value['status'] in ('idle', 'running', 'interrupted', 'complete')
                   and (value['phase'] is None or value['phase'] in PHASES325)
                   and (value['signal'] is None or
                        type(value['signal']) is int and value['signal'] > 0), 'session_schema')

    def _read(self, name):
        path = self.root / name
        require325(path.is_file() and not path.is_symlink() and
                   path.stat().st_size <= self.metadata_limit, 'metadata_file')
        try:
            # Duplicate keys must not be silently overwritten by JSON decoding.
            def object_pairs(pairs):
                result = {}
                for key, value in pairs:
                    require325(key not in result, 'duplicate_json_key')
                    result[key] = value
                return result
            return json.loads(path.read_bytes(), object_pairs_hook=object_pairs,
                              parse_constant=lambda x: require325(False, 'nonfinite_json'))
        except (UnicodeError, json.JSONDecodeError) as error:
            raise ResumeStateError325('resume325.metadata_json') from error

    def _atomic(self, name, value):
        data = canonical325(value)
        require325(len(data) <= self.metadata_limit, 'metadata_bound')
        temp = self.root / '.atomic.tmp'
        require325(not temp.exists() and not temp.is_symlink(), 'atomic_temp_exists')
        descriptor = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        try:
            with os.fdopen(descriptor, 'wb') as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temp, self.root / name)
            _fsync_dir(self.root)
        except BaseException:
            # Keep a partial temporary file after an interruption. It has no
            # authority and will be discarded only at an explicit resume.
            raise

    def _save_reservations(self):
        require325(sum(self._reservations.values()) <= self.payload_cap, 'reservation_cap')
        self._atomic('reservations.json', {'schema': 1, 'binding_sha256': self.binding_sha256,
                                        'allocations': self._reservations})

    def _save_session(self):
        self._check_session()
        self._atomic('session.json', self._session)

    def _inventory(self):
        total = 0
        for base, directories, files in os.walk(self.root, followlinks=False):
            base = Path(base)
            for name in directories + files:
                path = base / name
                metadata = path.lstat()
                require325(not stat.S_ISLNK(metadata.st_mode), 'symlink_in_store')
                require325(stat.S_ISDIR(metadata.st_mode) or stat.S_ISREG(metadata.st_mode),
                           'nonregular_store_entry')
                if stat.S_ISREG(metadata.st_mode):
                    require325(metadata.st_nlink == 1, 'hardlink_in_store')
                    total += metadata.st_size
                    if name in _ROOT_FILES or name == 'certificate.json':
                        require325(metadata.st_size <= self.metadata_limit, 'metadata_bound')
                rel = path.relative_to(self.root)
                if len(rel.parts) == 1:
                    require325((name in _ROOT_FILES and path.is_file()) or
                               (name in ('phases', 'partial') and path.is_dir()),
                               'unknown_root_entry')
                elif len(rel.parts) == 2:
                    require325(rel.parts[0] in ('phases', 'partial') and
                               name in PHASES325 and path.is_dir(), 'unknown_phase_path')
                else:
                    require325(rel.parts[0] in ('phases', 'partial') and
                               rel.parts[1] in PHASES325 and
                               all(_NAME.fullmatch(x) for x in rel.parts[2:]),
                               'unknown_artifact_path')
        require325(total <= self.cap_bytes, 'actual_root_cap')
        return total

    def _phase_files(self, directory):
        return {str(p.relative_to(directory)): p for p in directory.rglob('*') if p.is_file()}

    def validate_all(self):
        """Rehash *every* completed artifact, including complete SQLite bytes."""
        require325(self._lock is not None, 'closed_store')
        self._inventory()
        self.certificates = {}
        partial_names = {p.name for p in (self.root / 'partial').iterdir()}
        complete_names = {p.name for p in (self.root / 'phases').iterdir()}
        require325(len(partial_names) <= 1 and not (partial_names & complete_names),
                   'ambiguous_partial_phase')
        for phase, dependencies in PHASES325.items():
            directory = self.root / 'phases' / phase
            if not directory.exists():
                continue
            path = directory / 'certificate.json'
            cert = self._read(str(path.relative_to(self.root)))
            _fields(cert, ('schema', 'phase', 'binding_sha256', 'dependencies',
                           'artifacts', 'census', 'notes'), 'certificate')
            require325(type(cert['schema']) is int and cert['schema'] == 1 and cert['phase'] == phase and
                       cert['binding_sha256'] == self.binding_sha256, 'certificate_binding')
            require325(all(d in self.certificates for d in dependencies), 'missing_dependency')
            expected = {d: self.certificates[d]['certificate_sha256'] for d in dependencies}
            require325(cert['dependencies'] == expected, 'dependency_certificate_changed')
            self._check_evidence(cert['census'], cert['notes'])
            require325(type(cert['artifacts']) is dict and bool(cert['artifacts']),
                       'artifact_inventory')
            files = self._phase_files(directory)
            paths = []
            payload = 0
            for name, artifact in cert['artifacts'].items():
                require325(type(name) is str and _NAME.fullmatch(name), 'artifact_name')
                _fields(artifact, ('path', 'size', 'sha256'), 'artifact')
                rel = _relative(artifact['path'])
                paths.append(str(rel))
                item = directory / rel
                require325(type(artifact['size']) is int and artifact['size'] >= 0 and
                           type(artifact['sha256']) is str and _SHA.fullmatch(artifact['sha256']),
                           'artifact_schema')
                require325(item.is_file() and item.stat().st_size == artifact['size'] and
                           sha256325(item) == artifact['sha256'], 'artifact_bytes_changed')
                payload += artifact['size']
            require325(len(paths) == len(set(paths)) and
                       set(files) == set(paths) | {'certificate.json'}, 'artifact_census')
            require325(payload <= self._reservations.get(phase, -1), 'unreserved_artifact_bytes')
            self.certificates[phase] = {'certificate': cert,
                                        'certificate_sha256': sha256325(path)}
        for partial in (self.root / 'partial').iterdir():
            used = sum(p.stat().st_size for n, p in self._phase_files(partial).items()
                       if n != 'certificate.json')
            require325(used <= self._reservations.get(partial.name, 0),
                       'partial_allocation_exceeded')
        require325(sum(self._reservations.values()) <= self.payload_cap, 'reservation_cap')
        return copy.deepcopy(self.certificates)

    @staticmethod
    def _check_evidence(census, notes):
        require325(type(census) is dict and bool(census) and
                   type(notes) is list and all(type(n) is str for n in notes), 'phase_evidence')
        canonical325({'census': census, 'notes': notes})

    def completed(self, phase):
        require325(phase in PHASES325, 'unknown_phase')
        return phase in self.certificates

    def artifact(self, phase, name):
        require325(self.completed(phase), 'phase_not_complete')
        artifacts = self.certificates[phase]['certificate']['artifacts']
        require325(name in artifacts, 'unknown_artifact')
        return self.root / 'phases' / phase / artifacts[name]['path']

    def start(self, phase, *, discard_partial=False):
        self._writable()
        require325(phase in PHASES325 and self._lock is not None, 'unknown_phase_or_closed')
        require325(not self.completed(phase), 'phase_already_complete')
        require325(all(self.completed(d) for d in PHASES325[phase]), 'missing_dependency')
        self._inventory()
        partials = list((self.root / 'partial').iterdir())
        require325(all(p.name == phase for p in partials), 'other_partial_phase')
        stage = self.root / 'partial' / phase
        if stage.exists():
            require325(discard_partial is True and self.explicit_resume is True,
                       'explicit_partial_discard_required')
            # Only known staging files are removed, under the locked, bound,
            # independently protected workspace. No rmtree/caller path deletion.
            for base, directories, files in os.walk(stage, topdown=False, followlinks=False):
                for name in files:
                    (Path(base) / name).unlink()
                for name in directories:
                    (Path(base) / name).rmdir()
            stage.rmdir()
            _fsync_dir(stage.parent)
            self._reservations.pop(phase, None)
            self._save_reservations()
        stage.mkdir(mode=0o700)
        _fsync_dir(stage.parent)
        self._session = {'schema': 1, 'binding_sha256': self.binding_sha256,
                         'status': 'running', 'phase': phase, 'signal': None}
        self._save_session()
        return stage

    def reserve(self, phase, bytes_total):
        """Reserve a total payload allowance *before* external phase writes."""
        self._writable()
        require325(type(bytes_total) is int and bytes_total >= 0, 'reservation_size')
        require325(self._session['status'] == 'running' and self._session['phase'] == phase
                   and (self.root / 'partial' / phase).is_dir(), 'phase_not_running')
        self._inventory()
        used = sum(p.stat().st_size for n, p in self._phase_files(
            self.root / 'partial' / phase).items() if n != 'certificate.json')
        require325(used <= bytes_total, 'reservation_below_existing_bytes')
        updated = dict(self._reservations, **{phase: bytes_total})
        require325(sum(updated.values()) <= self.payload_cap, 'reservation_cap')
        self._reservations = updated
        self._save_reservations()
        return bytes_total

    def commit(self, phase, artifacts, *, census, notes=None):
        self._writable()
        require325(self._session['status'] == 'running' and self._session['phase'] == phase,
                   'phase_not_running')
        require325(type(artifacts) is dict and bool(artifacts), 'artifact_inventory')
        self._check_evidence(census, [] if notes is None else notes)
        self._inventory()
        stage = self.root / 'partial' / phase
        paths = {name: _relative(value) for name, value in artifacts.items()}
        require325(all(type(n) is str and _NAME.fullmatch(n) for n in paths), 'artifact_name')
        require325(len(set(paths.values())) == len(paths), 'duplicate_artifact_path')
        require325(set(self._phase_files(stage)) == {str(p) for p in paths.values()},
                   'artifact_census')
        proof = {}
        for name, rel in paths.items():
            path = stage / rel
            # The caller has closed the writer; flush and bind the complete bytes.
            with path.open('rb') as stream:
                os.fsync(stream.fileno())
            proof[name] = {'path': str(rel), 'size': path.stat().st_size,
                           'sha256': sha256325(path)}
        require325(sum(p['size'] for p in proof.values()) <=
                   self._reservations.get(phase, -1), 'unreserved_artifact_bytes')
        dependencies = {d: self.certificates[d]['certificate_sha256'] for d in PHASES325[phase]}
        certificate = {'schema': 1, 'phase': phase, 'binding_sha256': self.binding_sha256,
                       'dependencies': dependencies, 'artifacts': proof,
                       'census': census, 'notes': [] if notes is None else notes}
        self._atomic(str((stage / 'certificate.json').relative_to(self.root)), certificate)
        for base, _, _ in os.walk(stage, topdown=False):
            _fsync_dir(base)
        destination = self.root / 'phases' / phase
        require325(not destination.exists(), 'phase_already_complete')
        os.rename(stage, destination)
        _fsync_dir(destination.parent)
        _fsync_dir(stage.parent)
        self._reservations[phase] = sum(p['size'] for p in proof.values())
        self._save_reservations()
        self._session = {'schema': 1, 'binding_sha256': self.binding_sha256,
                         'status': 'complete', 'phase': phase, 'signal': None}
        self._save_session()
        self.validate_all()
        return copy.deepcopy(self.certificates[phase])

    def mark_interrupted(self, signum):
        self._writable()
        require325(isinstance(signum, int) and not isinstance(signum, bool) and signum > 0,
                   'signal_number')
        require325(not (self.root / 'semantic-failure.json').exists(), 'semantic_failed_latch')
        self._session = dict(self._session, status='interrupted', signal=int(signum))
        self._save_session()

    def mark_semantic_failed(self, detail):
        self._writable()
        require325(type(detail) is str and bool(detail), 'failure_detail')
        require325(not (self.root / 'semantic-failure.json').exists(), 'semantic_failed_latch')
        value = {'schema': 1, 'status': 'semantic_failed',
                 'binding_sha256': self.binding_sha256,
                 'phase': self._session['phase'], 'detail': detail}
        data = canonical325(value)
        require325(len(data) <= self.metadata_limit, 'metadata_bound')
        # Create the immutable failure latch *before* writing its body. A kill
        # while recording failure leaves a malformed latch that still blocks
        # resume, rather than an ordinary atomic temp that might be discarded.
        path = self.root / 'semantic-failure.json'
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        _fsync_dir(self.root)
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        _fsync_dir(self.root)
