#!/usr/bin/env python3
"""Fresh, read-only receipt for completed Stage325 job55576.

Run with python3 -B after terminal collection. This imports the separately
sealed native adapter solely for fresh input/source/checker binding. It never
opens PhaseStore for writing and never calls replay/provider/solver routines.
The binding method necessarily re-reads its two saved rank0 family inventories.
Only stdout is written; redirect it outside every protected input/state tree.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time

sys.dont_write_bytecode = True
JOB = 55576
EXPECTED_PACKAGE = '1a47b001cad28d82de231623f3c85e7b271e744369f9549d8a7b6b30d56dd196'
KIT = Path('/1/home/kevlahan/rk325-native-qualification-kit-v5p1')
STATE = Path.home() / 'rk325-native-v5p1-latest.json'
PHASES = {'cache': (), 'providers-C0': ('cache',), 'providers-L256': ('cache',),
          'pair': ('cache', 'providers-C0', 'providers-L256'),
          'work': ('pair',), 'preservation': ('work',)}
LIMIT = 1024**2


def require(ok, reason):
    if not ok:
        raise ValueError('postterminal55576.' + reason)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                      allow_nan=False).encode('ascii') + b'\n'


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def decode(data):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate_json_field:' + key[:128])
            result[key] = value
        return result
    return json.loads(data, object_pairs_hook=pairs,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def info(path):
    value = Path(path).stat()
    return dict(device=value.st_dev, inode=value.st_ino, bytes=value.st_size,
                mtime_ns=value.st_mtime_ns, ctime_ns=value.st_ctime_ns,
                links=value.st_nlink, owner=value.st_uid)


def regular(path):
    path = Path(path)
    require(path.is_absolute() and '..' not in path.parts, 'absolute_clean_path')
    for component in (path, *path.parents):
        require(not component.is_symlink(), 'symlink:' + str(component))
    require(stat.S_ISREG(path.lstat().st_mode), 'regular_file:' + str(path))
    return path


def read_json(path):
    path = regular(path)
    before = info(path)
    require(before['bytes'] <= LIMIT, 'json_bound:' + str(path))
    raw = path.read_bytes()
    require(info(path) == before and len(raw) == before['bytes'], 'json_changed:' + str(path))
    return decode(raw), hashlib.sha256(raw).hexdigest()


def hash_stable(path):
    path = regular(path)
    before = info(path)
    require(before['links'] == 1, 'database_or_artifact_hardlink:' + str(path))
    count = 0
    value = hashlib.sha256()
    with path.open('rb') as source:
        while chunk := source.read(4 * 1024**2):
            count += len(chunk)
            value.update(chunk)
    after = info(path)
    require(after == before and count == before['bytes'], 'file_changed_during_hash:' + str(path))
    return dict(path=str(path), bytes=count, sha256=value.hexdigest(),
                stat_before=before, stat_after=after)


def accounting():
    command = ['sacct', '-n', '-P', '-X', '-j', str(JOB), '-o',
               'JobIDRaw,State,ExitCode,ElapsedRaw,Start,End,NodeList']
    process = subprocess.run(command, capture_output=True, text=True, timeout=30, check=False)
    require(process.returncode == 0 and len(process.stdout) + len(process.stderr) <= 65536,
            'scheduler_failed_or_oversized')
    rows = [line.rstrip('|').split('|') for line in process.stdout.splitlines() if line.strip()]
    require(len(rows) == 1 and len(rows[0]) == 7 and rows[0][0] == str(JOB),
            'scheduler_missing_or_ambiguous_job')
    row = rows[0]
    require(row[1] == 'COMPLETED' and row[2] == '0:0' and row[6] == 'bb13',
            'job_not_cleanly_completed')
    require(row[3].isdigit() and 0 < int(row[3]) <= 1800 and
            row[4] not in ('Unknown', '', 'None') and row[5] not in ('Unknown', '', 'None'),
            'scheduler_elapsed_or_terminal_timestamp')
    return dict(command=command, returncode=0, row=row, stdout=process.stdout,
                stderr=process.stderr, terminal=True, clean_completion=True)


def safe_relative(value):
    require(type(value) is str, 'relative_path_type')
    path = Path(value)
    require(not path.is_absolute() and '..' not in path.parts and value == str(path),
            'unsafe_relative_path')
    return path


def main():
    started = time.monotonic()
    require(len(sys.argv) == 1, 'fixed_job_no_arguments')
    before_account = accounting()
    require(not STATE.is_symlink(), 'submission_state_symlink')
    state, state_sha = read_json(STATE.resolve(strict=True))
    require(state['job'] == JOB and state['source_job'] == 55569 and
            state['package_manifest_sha256'] == EXPECTED_PACKAGE and
            state['automatic_retry'] is False and
            state['builds'] == state['solver_runs'] == state['performance_runs'] == 0,
            'submission_identity_or_scope')
    kit = KIT.resolve(strict=True)
    require(Path(state['kit']).resolve(strict=True) == kit, 'submitted_kit_differs')
    package_path = regular(kit / 'native-files-sha256.json')
    require(hashlib.sha256(package_path.read_bytes()).hexdigest() == EXPECTED_PACKAGE,
            'sealed_package_manifest_differs')
    # No package code is imported until its pinned manifest has been checked.
    manifest, _ = read_json(package_path)
    actual = {}
    for base, directories, files in os.walk(kit, followlinks=False):
        for name in directories + files:
            path = Path(base) / name
            metadata = path.lstat()
            require(not stat.S_ISLNK(metadata.st_mode) and
                    (stat.S_ISDIR(metadata.st_mode) or stat.S_ISREG(metadata.st_mode)),
                    'package_entry_type')
            if stat.S_ISREG(metadata.st_mode):
                actual[str(path.relative_to(kit))] = path
    require(set(actual) == set(manifest) | {'native-files-sha256.json'}, 'package_inventory')
    for name, expected in manifest.items():
        require(bool(re.fullmatch('[0-9a-f]{64}', expected)) and
                hash_stable(actual[str(safe_relative(name))])['sha256'] == expected,
                'package_payload_differs:' + name)
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
        os.environ[name] = '1'
    sys.path.insert(0, str(kit))
    from native_runner325 import package_check325
    from resume_native325 import NativeBackend325
    require(package_check325(kit) == EXPECTED_PACKAGE, 'native_package_checker_differs')
    lexical_out = Path(state['out'])
    require(lexical_out.is_absolute() and '..' not in lexical_out.parts, 'output_root_path')
    canonical_out = lexical_out.resolve(strict=True)
    require(canonical_out.is_dir() and os.path.samefile(lexical_out, canonical_out),
            'output_root_alias')
    left_stat, right_stat = lexical_out.stat(), canonical_out.stat()
    output_identity = dict(lexical_path=str(lexical_out), canonical_path=str(canonical_out),
        samefile=os.path.samefile(lexical_out, canonical_out),
        lexical_stat=dict(device=left_stat.st_dev, inode=left_stat.st_ino),
        canonical_stat=dict(device=right_stat.st_dev, inode=right_stat.st_ino))
    state_root = canonical_out / 'state'
    require(not (state_root / 'semantic-failure.json').exists(), 'semantic_failure_latch')
    require((state_root / 'partial').is_dir() and not any((state_root / 'partial').iterdir()),
            'uncommitted_partial_state')
    identity, _ = read_json(state_root / 'identity.json')
    wrapped_binding = identity['binding']
    require(identity['binding_sha256'] == digest(wrapped_binding) and
            identity['root'] == str(state_root), 'state_identity_binding')
    completion, _ = read_json(canonical_out / 'worker-complete.json')
    attempt, _ = read_json(canonical_out / 'state.attempt.json')
    supervisor, _ = read_json(canonical_out / 'supervisor-report.json')
    require(completion['passed'] is True and attempt['status'] == 'passed' and
            attempt['native_gate_passed'] is True and supervisor['status'] == 'passed' and
            supervisor['gate_passed'] is True, 'complete_native_gate_required')
    certificate_hashes = {}
    certificates = {}
    for phase, dependencies in PHASES.items():
        certificate, certificate_sha = read_json(state_root / 'phases' / phase / 'certificate.json')
        require(set(certificate) == {'schema','phase','binding_sha256','dependencies','artifacts','census','notes'}
                and certificate['schema'] == 1 and certificate['phase'] == phase and
                certificate['binding_sha256'] == digest(wrapped_binding) and
                certificate['census']['passed'] is True,
                'closed_phase_certificate:' + phase)
        require(certificate['dependencies'] == {name: certificate_hashes[name] for name in dependencies},
                'phase_dependency:' + phase)
        certificate_hashes[phase] = certificate_sha
        certificates[phase] = certificate
    require(certificate_hashes == completion['phase_certificates'] == attempt['certificate_sha256'] ==
            supervisor['phase_certificates'], 'six_certificate_maps_differ')
    require(set(certificates['cache']['artifacts']) == {'left_db','right_db','receipt'},
            'cache_artifact_census')
    preservation_rel = safe_relative(certificates['preservation']['artifacts']['report']['path'])
    preservation, preservation_sha = read_json(state_root / 'phases' / 'preservation' / preservation_rel)
    require(preservation_sha == certificates['preservation']['artifacts']['report']['sha256'] and
            preservation['passed'] is True and preservation['full_input_hashes_matched'] == 1536,
            'preservation_artifact')
    # Fresh adapter construction uses only the validated sealed package. Its
    # binding() rehashes all1536 historical inputs, all original metadata/code,
    # binary/checkpoint records, and the checked-source reuse certificates.
    backend = NativeBackend325(kit / 'recovery_v4p1')
    fresh_binding = backend.binding()
    require(fresh_binding == preservation['binding'], 'fresh_raw_binding_differs')
    require(len(backend.anchors) == 1536 and backend._bound_metadata['passed'] is True and
            backend._bound_metadata['source_files'] == 258 and
            backend._bound_metadata['original_kit_files'] == 79 and
            backend._bound_metadata['metadata_files'] == 3889,
            'fresh_preservation_census')
    source_binding = backend.bind_sources()
    require(source_binding == preservation['source_binding'] and
            backend._bound_metadata == preservation['metadata'] and
            backend._checked_reuse == preservation['checked_reuse'], 'fresh_source_or_reuse_differs')
    databases = {}
    for name in ('left_db', 'right_db'):
        expected = certificates['cache']['artifacts'][name]
        require(set(expected) == {'path','size','sha256'}, 'database_certificate_schema')
        relative = safe_relative(expected['path'])
        path = state_root / 'phases' / 'cache' / relative
        require(path.parent.resolve(strict=True).is_relative_to(state_root / 'phases' / 'cache'),
                'database_outside_cache')
        for suffix in ('-wal', '-shm', '-journal'):
            require(not os.path.lexists(str(path) + suffix), 'database_transaction_sidecar')
        row = hash_stable(path)
        require(row['bytes'] == expected['size'] and row['sha256'] == expected['sha256'],
                'retained_database_bytes_differ:' + name)
        for suffix in ('-wal', '-shm', '-journal'):
            require(not os.path.lexists(str(path) + suffix), 'database_transaction_sidecar_after_hash')
        lexical_path = lexical_out / 'state' / 'phases' / 'cache' / relative
        require(os.path.samefile(path, lexical_path), 'database_alias_differs')
        row.update(canonical_path=str(path), lexical_path=str(lexical_path),
                   samefile=True, expected_bytes=expected['size'], expected_sha256=expected['sha256'])
        databases[name] = row
    require(package_check325(kit) == EXPECTED_PACKAGE, 'package_changed_during_post_terminal_hashing')
    after_account = accounting()
    require(after_account['row'] == before_account['row'], 'scheduler_identity_changed')
    require(canonical_out.stat().st_dev == right_stat.st_dev and
            canonical_out.stat().st_ino == right_stat.st_ino, 'output_root_changed')
    return dict(passed=True, job=JOB, package_manifest_sha256=EXPECTED_PACKAGE,
        scheduler=before_account, scheduler_after=after_account, state_file_sha256=state_sha,
        output_identity=output_identity, fresh_binding=fresh_binding,
        metadata=backend._bound_metadata, source_binding=source_binding,
        checked_reuse=backend._checked_reuse, full_input_hashes_matched=1536,
        phase_certificate_sha256=certificate_hashes, databases=databases,
        binding_scope='Fresh original metadata/source/checker/config plus1536 input hashes and2 intrinsic saved rank0 family-inventory parses; no full semantic/provider comparison or solver rerun.',
        elapsed_seconds=time.monotonic()-started, performance_qualified=False,
        native_qualification=False, independent_review_required=True,
        mutations='None: no PhaseStore lock/session/reserve/commit APIs; stdout only.')


if __name__ == '__main__':
    try:
        result = main()
    except BaseException as error:
        result = dict(passed=False, job=JOB, package_manifest_sha256=EXPECTED_PACKAGE,
                      error_type=type(error).__name__, error=str(error)[:4096],
                      native_qualification=False, performance_qualified=False)
    encoded = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n'
    require(len(encoded.encode()) <= LIMIT, 'final_stdout_report_bound')
    sys.stdout.write(encoded)
    raise SystemExit(0 if result['passed'] else 1)
