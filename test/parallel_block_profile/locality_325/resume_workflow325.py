"""Offline orchestration of closed Stage325 verification phases.

No scheduler, compiler or solver entry point is provided. A native invocation
requires a separately authorized time budget. Fixture completion never grants
native qualification. PhaseStore certificates, not this progress report, control
reuse. All callbacks write only inside their reserved private staging directory.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time

from resume_state325 import PhaseStore325, canonical325, require325, sha256325
import stream_compare325 as stream
import verify325 as q

MiB = 1024**2
REPORT_LIMIT325 = MiB
TABLES325 = ('route', 'domain', 'context', 'feature')
NATIVE_ROWS325 = {'route': 65024, 'domain': 256, 'context': 66048, 'feature': 508840}


class ReplayInterrupted325(Exception):
    def __init__(self, signum):
        self.signum = signum
        super().__init__(signal.Signals(signum).name)


class SemanticFailure325(ValueError):
    pass


def utc325():
    return datetime.now(timezone.utc).isoformat()


def save325(path, value, *, limit=REPORT_LIMIT325):
    data = canonical325(value)
    require325(len(data) <= limit, 'report_bound')
    path = Path(path)
    require325(not path.is_symlink(), 'report_symlink')
    with path.open('xb') as output:
        output.write(data)
        output.flush()
        os.fsync(output.fileno())


def _sha_value(value):
    return hashlib.sha256(canonical325(value)).hexdigest()


def _load(path):
    # Artifacts are already fully byte-bound by PhaseStore. Duplicate fields
    # are nevertheless rejected, rather than overwritten by JSON decoding.
    def pairs(items):
        result = {}
        for key, value in items:
            require325(key not in result, 'duplicate_report_field')
            result[key] = value
        return result
    return json.loads(Path(path).read_bytes(), object_pairs_hook=pairs)


def inspect_cache325(path, arm, rows, expected_ranks):
    """Exact completeness before certification; pair equality alone is insufficient."""
    require325(set(rows) == set(TABLES325) and all(type(x) is int and x > 0 for x in rows.values()),
               'expected_cache_census')
    with stream.readonly_database325(path) as (connection, identity):
        require325(connection.execute('pragma integrity_check').fetchall() == [('ok',)],
                   'cache_integrity')
        identity['schema_sha256'] = stream._schema325(connection, TABLES325)
        stream._arm325(connection, TABLES325, arm)
        actual = {table: connection.execute('select count(*) from ' + table).fetchone()[0]
                  for table in TABLES325}
        if actual != rows:
            raise SemanticFailure325('cache exact table census: ' + str(actual) + ' expected ' + str(rows))
        identity['domain_rank_census_sha256'] = stream._ranks325(connection, 'domain', arm, expected_ranks)
    return dict(arm=arm, rows=actual, identity=identity, sha256=sha256325(path))


def _normalize_paths(value, stage):
    """Keep reports valid after the staging directory's atomic rename."""
    if isinstance(value, dict):
        return {key: _normalize_paths(item, stage) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_paths(item, stage) for item in value]
    if isinstance(value, str) and value.startswith(str(stage) + '/'):
        return 'phase:/' + value[len(str(stage)) + 1:]
    return value


@contextmanager
def _signals325(seconds, report):
    require325(type(seconds) in (int, float) and 0 < seconds < 86400, 'explicit_time_budget')
    require325(signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0), 'existing_timer')
    originals = {}
    def stop(signum, frame):
        # Python delivers this in the main thread. Persist outside the handler;
        # no filesystem/SQLite calls run from interrupted writer code.
        report['signal'] = dict(number=signum, name=signal.Signals(signum).name,
                                utc=utc325(), monotonic=time.monotonic(), phase=report['phase'])
        raise ReplayInterrupted325(signum)
    try:
        for signum in (signal.SIGTERM, signal.SIGALRM, signal.SIGINT):
            originals[signum] = signal.signal(signum, stop)
        signal.setitimer(signal.ITIMER_REAL, seconds)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        for signum, handler in originals.items():
            signal.signal(signum, handler)


@contextmanager
def _held_store325(store, record_failure):
    with store:
        try:
            yield store
        except BaseException as exc:
            try:
                if isinstance(exc, ReplayInterrupted325):
                    store.mark_interrupted(exc.signum)
                elif isinstance(exc, (SemanticFailure325, q.v.AuditError)) and not (store.root/'semantic-failure.json').exists():
                    store.mark_semantic_failed(str(exc)[:4096])
            finally:
                record_failure(exc)
            raise


def run325(backend, root, *, max_seconds, explicit_resume=False, discard_partial=False,
           cap_bytes=4 * 1024**3, phase_limits=None):
    """Verify existing captures; return a bounded report or raise on interruption/failure.

    backend.binding() MUST freshly hash original inputs, code and configuration.
    NativeBackend325 supplies the production contract. Offline fixtures explicitly
    set native=False and reduced expected census; they cannot qualify native data.
    The current attempt report and its atomic temporary file reserve2MiB outside
    the store, within the total logical-byte cap. Completed payloads remain owned
    by the store. No cleanup occurs after a failure or signal.
    """
    root = Path(root)
    require325(root.is_absolute() and '..' not in root.parts, 'workflow_absolute_root')
    for part in (root, *root.parents):
        require325(not part.is_symlink(), 'workflow_symlink')
    protected = tuple(Path(x) for x in backend.protected_roots)
    report_path = root.with_name(root.name + '.attempt.json')
    temp = root.with_name(root.name + '.attempt.tmp')
    for item in (report_path, temp):
        require325(not item.is_symlink() and (not item.exists() or item.is_file()), 'attempt_file_type')
        require325(not item.exists() or (item.stat().st_nlink == 1 and item.stat().st_size <= REPORT_LIMIT325),
                   'attempt_file_bound_or_link')
        require325(all(not item.is_relative_to(p) and not p.is_relative_to(item) for p in protected),
                   'protected_report_overlap')
    limits = dict({'providers-C0': 2*MiB, 'providers-L256': 2*MiB, 'pair': 2*MiB,
                   'work': 192*MiB, 'preservation': 32*MiB}, **(phase_limits or {}))
    require325(set(limits) <= {'cache', 'providers-C0', 'providers-L256', 'pair', 'work', 'preservation'}
               and all(type(x) is int and x > 0 for x in limits.values()), 'phase_limits')
    ranks, steps, families = backend.expected_ranks, backend.steps, backend.expected_families
    arms = tuple(backend.arm_names)
    require325(arms == ('qual-C0', 'qual-L256'), 'workflow_arm_identity')
    expected_rows = backend.expected_table_rows
    if backend.native:
        from resume_native325 import NativeBackend325
        require325(type(backend) is NativeBackend325, 'native_backend_required')
        require325((ranks, steps, families) == (256, 4, 254) and expected_rows == NATIVE_ROWS325,
                   'native_exact_contract')
        require325(backend.expected_family_keys is None or len(backend.expected_family_keys) == 254,
                   'native_family_contract')
    started = time.monotonic()
    report = dict(schema=1, status='running', phase='fresh_input_binding', started_utc=utc325(),
                  started_monotonic=started, max_seconds=max_seconds, logical_byte_cap=cap_bytes,
                  native=backend.native, native_qualification=False, performance_qualified=False,
                  builds=0, solver_runs=0, performance_runs=0, signal=None, phases=[],
                  child_exit_statuses=[], retention='not_opened')
    owned = False
    store = None
    def persist():
        # Fixed paths are owned by this one locked store. A killed atomic report
        # writer leaves a recognized temporary file; it has no proof authority.
        report['elapsed_seconds'] = time.monotonic() - started
        statuses = getattr(backend, 'child_exit_statuses', [])
        require325(type(statuses) is list and len(statuses) <= 64, 'child_status_bound')
        report['child_exit_statuses'] = statuses
        data = canonical325(report)
        require325(len(data) <= REPORT_LIMIT325, 'attempt_report_bound')
        require325(not temp.is_symlink() and not report_path.is_symlink(), 'attempt_symlink')
        with temp.open('wb') as output:
            output.write(data); output.flush(); os.fsync(output.fileno())
        os.replace(temp, report_path)
        descriptor = os.open(root.parent, os.O_RDONLY | os.O_DIRECTORY)
        try: os.fsync(descriptor)
        finally: os.close(descriptor)

    def execute(phase, callback):
        report['phase'] = phase
        timing = dict(phase=phase, started_utc=utc325(), started_monotonic=time.monotonic(), reused=False)
        report['phases'].append(timing)
        persist()

        if store.completed(phase):
            timing.update(reused=True, elapsed_seconds=time.monotonic()-timing['started_monotonic'])
            persist()
            return
        stage = store.start(phase, discard_partial=discard_partial)
        if phase == 'cache' and phase not in limits:
            allowance = store.available_bytes() - sum(limits.values())
        else:
            allowance = limits[phase]
        require325(allowance > 0, 'no_phase_space')
        store.reserve(phase, allowance)
        timing['allowance_bytes'] = allowance
        callback_start = time.monotonic()
        result = callback(stage, allowance)
        timing['callback_seconds'] = time.monotonic() - callback_start
        require325(type(result) is dict and set(result) == {'artifacts', 'census'} and
                   result['census'].get('passed') is True, 'incomplete_phase_result')
        commit_start = time.monotonic()
        store.commit(phase, result['artifacts'], census=result['census'])
        timing['commit_and_full_rehash_seconds'] = time.monotonic() - commit_start
        timing['elapsed_seconds'] = time.monotonic() - timing['started_monotonic']
        timing['certificate_sha256'] = store.certificates[phase]['certificate_sha256']
        report['retention'] = 'closed_certificates_and_uncommitted_partial_retained'
        persist()

    def record_failure(exc):
        report.update(status='interrupted' if isinstance(exc, ReplayInterrupted325) else 'failed',
                      error=dict(type=type(exc).__name__, message=str(exc)[:4096]), completed_utc=utc325())
        if report['phases']:
            current = report['phases'][-1]
            current.update(status=report['status'], elapsed_seconds=time.monotonic()-current['started_monotonic'])
        persist()

    try:
        with _signals325(max_seconds, report):
            binding = backend.binding()
            report['fresh_binding_seconds'] = time.monotonic() - started
            require325(backend.expected_family_keys is not None and len(backend.expected_family_keys) == families,
                       'bound_expected_family_contract')
            contract = dict(ranks=ranks, steps=steps, families=families, arms=arms,
                            rows=expected_rows, native=backend.native, limits=limits,
                            family_keys=backend.expected_family_keys)
            binding = dict(binding, config=_sha_value(dict(backend=binding['config'], contract=contract)))
            # Bind the orchestrator and both new helpers even if a backend omits
            # them accidentally. NativeBackend additionally binds all old gates.
            binding['checker'] = _sha_value(dict(backend=binding['checker'], files={
                name: sha256325(Path(__file__).with_name(name)) for name in
                ('resume_workflow325.py', 'resume_state325.py', 'stream_compare325.py')}))
            store = PhaseStore325(root, binding, protected_roots=protected,
                                  cap_bytes=cap_bytes-2*REPORT_LIMIT325, explicit_resume=explicit_resume)
            report['phase'] = 'validate_retained_certificates'
            validation_started = time.monotonic()
            with _held_store325(store, record_failure):
                owned = True
                report['retained_validation_seconds'] = time.monotonic() - validation_started
                report['binding'] = binding
                if report_path.exists():
                    prior = _load(report_path)
                    require325(prior.get('schema') == 1 and prior.get('binding') == binding,
                               'attempt_history_binding')
                    history = prior.get('attempt_history', [])
                    require325(type(history) is list and len(history) < 63, 'attempt_history_bound')
                    summary = {key: prior.get(key) for key in ('started_utc','completed_utc','status','phase',
                        'signal','error','elapsed_seconds','phases','child_exit_statuses')}
                    report['attempt_history'] = history + [summary]
                else:
                    report['attempt_history'] = []
                report['retention'] = 'closed_certificates_and_uncommitted_partial_retained'
                backend.bind_sources()
                def cache(stage, allowance):
                    result = backend.cache(stage, allowance)
                    if (type(result) is not dict or set(result) != {'artifacts','census'} or
                        result['census'].get('passed') is not True or
                        result['census'].get('streams') != 2*ranks or
                        result['census'].get('capture_reads') != 2*ranks):
                        raise SemanticFailure325('cache stream/actual-read closure incomplete')
                    census = []
                    for arm, key in zip(arms, ('left_db', 'right_db')):
                        require325(key in result['artifacts'], 'cache_database_artifact')
                        rel = Path(result['artifacts'][key])
                        require325(not rel.is_absolute() and '..' not in rel.parts, 'cache_artifact_path')
                        census.append(inspect_cache325(stage/rel, arm, expected_rows, ranks))
                    result['census'] = _normalize_paths(dict(result['census'], exact_databases=census), stage)
                    return result
                execute('cache', cache)
                dbs = [store.artifact('cache', key) for key in ('left_db', 'right_db')]
                for arm, phase, db in zip(arms, ('providers-C0', 'providers-L256'), dbs):
                    def provider(stage, allowance, arm=arm, db=db):
                        def progress(stage_name, progress_arm, completed, total, detail=None, save=False):
                            report['progress_only'] = dict(stage=stage_name, arm=progress_arm,
                                completed=completed, total=total, detail=detail)
                            if save: persist()
                        certificate = stream.certify_providers325(db, arm, expected_ranks=ranks,
                            expected_families=families, expected_family_keys=backend.expected_family_keys,
                            progress=progress)
                        # No feature export intervenes between proof completion
                        # and durable publication of this compact result.
                        save325(stage/'providers.json', certificate, limit=allowance)
                        return dict(artifacts={'proof':'providers.json'}, census=dict(passed=True,
                            families=families, database_sha256=sha256325(db)))
                    execute(phase, provider)
                def pair(stage, allowance):
                    providers = [_load(store.artifact(phase, 'proof')) for phase in ('providers-C0','providers-L256')]
                    provider_pair = stream.compare_provider_certificates325(*providers,
                        expected_ranks=ranks, expected_families=families)
                    physical = stream.compare_databases325(*dbs, expected_ranks=ranks, steps=steps)
                    result = dict(passed=provider_pair['passed'] and physical['passed'],
                                  providers=provider_pair, physical=physical)
                    if not result['passed']:
                        # Publish the fail-closed latch before any diagnostic
                        # write can be interrupted. A partial report is not proof.
                        store.mark_semantic_failed('provider or physical cross-arm difference')
                        save325(stage/'pair.json', result, limit=allowance)
                        raise SemanticFailure325('provider or physical cross-arm difference; bounded pair.json retained')
                    save325(stage/'pair.json', result, limit=allowance)
                    return dict(artifacts={'proof':'pair.json'}, census=dict(passed=True,
                        feature_rows=physical['features']['sides'], bridges=physical['inverse_bridge']['proved']))
                execute('pair', pair)
                execute('work', backend.work)
                execute('preservation', backend.preservation)
                report.update(status='passed', phase='complete',
                              native_gate_passed=bool(backend.native), native_qualification=False,
                              independent_review_required=bool(backend.native), completed_utc=utc325(),
                              retention='all_six_closed_phase_certificates_retained',
                              certificate_sha256={name: row['certificate_sha256']
                                                  for name,row in store.certificates.items()})
                persist()
    except BaseException as exc:
        report.update(status='interrupted' if isinstance(exc, ReplayInterrupted325) else 'failed',
                      error=dict(type=type(exc).__name__, message=str(exc)[:4096]), completed_utc=utc325())
        if not owned:
            # Before the store lock is held, do not overwrite another writer's
            # attempt report. A bounded stderr event still identifies a signal
            # during input binding or certificate rehash, with its exact phase.
            event = {key: report[key] for key in ('status','phase','signal','error','completed_utc')}
            sys.stderr.write(canonical325(dict(preownership_termination=event)).decode('ascii'))
            sys.stderr.flush()
        raise
    return report
