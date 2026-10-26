"""Read-only, bounded Stage325 provider and cross-arm verification.

This module never exports feature dictionaries. Callers must bind the original
inverse sources before use and bind a completed database to its full SHA256 in
persistent phase state. File stat/schema identities here are additional guards,
not a substitute for that external full-file hash.
"""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
from urllib.parse import quote

import inverse_identity325 as inv
import verify325 as q

CONTRACT325 = 'stage325-stream-comparison-v1'
MAX_EXAMPLES325 = 32
COLUMNS325 = {
    'feature': ('arm', 'key', 'digest'),
    'context': ('arm', 'rank', 'kind', 'key', 'payload'),
    'route': ('arm', 'rank', 'kind', 'key', 'payload'),
    'domain': ('arm', 'rank', 'payload'),
}


def _require(ok, code, detail=None):
    if not ok:
        raise q.v.AuditError(code + (': ' + str(detail) if detail is not None else ''))


def _encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode() + b'\n'


def _digest(value):
    return hashlib.sha256(_encoded(value)).hexdigest()


class _Events325:
    def __init__(self):
        self.counts = {}
        self.examples = []
        self.hash = hashlib.sha256()

    def add(self, code, **detail):
        self.counts[code] = self.counts.get(code, 0) + 1
        event = dict(code=code, **detail)
        self.hash.update(_encoded(event))
        if len(self.examples) < MAX_EXAMPLES325:
            # Payloads and raw bridge contexts never enter bounded diagnostics.
            self.examples.append({key: (value[:512] if isinstance(value, str) else value)
                                  for key, value in event.items()})

    def report(self):
        return dict(count=sum(self.counts.values()), counts=self.counts,
                    sha256=self.hash.hexdigest(), examples=self.examples,
                    examples_limit=MAX_EXAMPLES325)


def _path325(database):
    _require('..' not in Path(database).parts, 'stream.database_parent_component')
    path = Path(os.path.abspath(os.fspath(database)))
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        info = current.lstat()
        _require(not stat.S_ISLNK(info.st_mode), 'stream.database_symlink', str(current))
        _require(stat.S_ISREG(info.st_mode) if current == path else stat.S_ISDIR(info.st_mode),
                 'stream.database_file_type', str(current))
    for suffix in ('-wal', '-shm', '-journal'):
        _require(not os.path.lexists(str(path) + suffix), 'stream.database_sidecar', suffix)
    return path


def _identity325(path):
    value = path.stat()
    return dict(device=value.st_dev, inode=value.st_ino, bytes=value.st_size,
                mtime_ns=value.st_mtime_ns, ctime_ns=value.st_ctime_ns)


@contextmanager
def readonly_database325(database):
    """Open without writes; reject unresolved transaction sidecars and symlinks."""
    path = _path325(database)
    before = _identity325(path)
    with path.open('rb') as source:
        header = source.read(20)
    _require(len(header) == 20 and header[:16] == b'SQLite format 3\x00', 'stream.database_header')
    _require(header[18:20] == b'\x01\x01', 'stream.database_wal_mode')
    # Even opening a WAL-mode file read-only can create -wal/-shm files.
    # Reject its header before opening, and immutable mode forbids such writes.
    uri = 'file:' + quote(str(path), safe='/') + '?mode=ro&immutable=1'
    connection = sqlite3.connect(uri, uri=True)
    try:
        connection.execute('pragma query_only=on')
        _require(connection.execute('pragma query_only').fetchone() == (1,), 'stream.query_only')
        _require(connection.execute('pragma journal_mode').fetchone()[0].lower() != 'wal',
                 'stream.database_wal_mode')
        connection.execute('begin')
        yield connection, dict(path=str(path), **before)
    finally:
        connection.close()
        _path325(path)
        _require(_identity325(path) == before, 'stream.database_changed', str(path))


def _schema325(connection, tables):
    schema = {}
    for table in tables:
        kind = connection.execute('select type from sqlite_master where name=?', (table,)).fetchall()
        _require(kind == [('table',)], 'stream.schema_table', table)
        rows = connection.execute('pragma table_info(' + table + ')').fetchall()
        _require(tuple(row[1] for row in rows) == COLUMNS325[table], 'stream.schema_columns', table)
        schema[table] = rows
    return _digest(schema)


def _arm325(connection, tables, arm):
    for table in tables:
        _require(connection.execute('select count(*) from ' + table + ' where arm is null or arm!=?',
                                    (arm,)).fetchone()[0] == 0, 'stream.unexpected_arm', table)


def _ranks325(connection, table, arm, expected_ranks):
    rows = connection.execute('select rank from ' + table + ' where arm=? order by rank', (arm,))
    digest = hashlib.sha256()
    for rank in range(expected_ranks):
        row = next(rows, None)
        _require(row == (rank,), 'stream.rank_census', dict(table=table, expected=rank, actual=row))
        digest.update(_encoded(row))
    _require(next(rows, None) is None, 'stream.extra_rank', table)
    return digest.hexdigest()


def _feature_rows325(connection, arm, events, side, stats):
    previous = None
    cursor = connection.execute('select key,digest from feature where arm=? order by key collate binary', (arm,))
    for key, digest in cursor:
        _require(isinstance(key, str) and isinstance(digest, str) and len(digest) == 64 and
                 all(c in '0123456789abcdef' for c in digest), 'stream.feature_record', side)
        stats['rows'] += 1
        stats['hash'].update(_encoded((key, digest)))
        if key == previous:
            events.add('duplicate_feature', side=side, key=key)
            continue
        previous = key
        stats['unique'] += 1
        yield key, digest
    if not stats['rows']:
        events.add('empty_features', side=side)


def _bridge_rows325(connection, arm, expected_ranks, steps, events, side, stats):
    previous = -1
    cursor = connection.execute("select rank,key,payload from context where arm=? and kind='inverse-bridge' "
                                'order by rank,cast(key as integer),key collate binary', (arm,))
    for rank, step_text, payload in cursor:
        _require(isinstance(rank, int) and 0 <= rank < expected_ranks and
                 isinstance(step_text, str) and step_text in tuple(str(x) for x in range(1, steps + 1)),
                 'stream.bridge_key', dict(side=side, rank=rank, step=step_text))
        step = int(step_text)
        index = rank * steps + step - 1
        _require(isinstance(payload, str), 'stream.bridge_payload', side)
        stats['rows'] += 1
        stats['hash'].update(_encoded((rank, step, payload)))
        if index == previous:
            events.add('duplicate_bridge', side=side, rank=rank, step=step)
            continue
        _require(index > previous, 'stream.bridge_order', side)
        for missing in range(previous + 1, index):
            events.add('missing_bridge', side=side, rank=missing // steps, step=missing % steps + 1)
        previous = index
        context = json.loads(payload)
        _require(type(context.get('step')) is int and context['step'] == step, 'stream.bridge_step_binding', dict(side=side, rank=rank, step=step))
        stats['unique'] += 1
        yield (rank, step), context
    for missing in range(previous + 1, expected_ranks * steps):
        events.add('missing_bridge', side=side, rank=missing // steps, step=missing % steps + 1)


def compare_databases325(left_db, right_db, *, expected_ranks=256, steps=4,
                         left_arm='qual-C0', right_arm='qual-L256'):
    """Compare sorted feature cursors and one strict physical bridge pair at a time.

    Semantic differences return passed=False with complete event counts/digests
    and at most32 examples. Malformed structure/path/source binding raises.
    Expected ranks/steps may be reduced only by offline fixtures; the native
    runner must hard-code256 and4.
    """
    _require(isinstance(expected_ranks, int) and expected_ranks > 0 and isinstance(steps, int) and steps > 0,
             'stream.configuration')
    inv.require_inverse_source_binding325()
    events = _Events325()
    feature_stats = [dict(rows=0, unique=0, hash=hashlib.sha256()) for _ in range(2)]
    bridge_stats = [dict(rows=0, unique=0, hash=hashlib.sha256()) for _ in range(2)]
    comparison_hash = hashlib.sha256()
    bridge_hash = hashlib.sha256()
    common = equal = bridge_pairs = bridge_proved = moved = 0
    with readonly_database325(left_db) as (left, lid), readonly_database325(right_db) as (right, rid):
        identities = []
        for connection, arm, identity in ((left, left_arm, lid), (right, right_arm, rid)):
            identity['schema_sha256'] = _schema325(connection, ('feature', 'context', 'domain'))
            _arm325(connection, ('feature', 'context', 'domain'), arm)
            identity['domain_rank_census_sha256'] = _ranks325(connection, 'domain', arm, expected_ranks)
            identities.append(identity)
        left_rows = _feature_rows325(left, left_arm, events, 'left', feature_stats[0])
        right_rows = _feature_rows325(right, right_arm, events, 'right', feature_stats[1])
        a, b = next(left_rows, None), next(right_rows, None)
        while a is not None or b is not None:
            if b is None or (a is not None and a[0] < b[0]):
                events.add('missing_feature_right', key=a[0])
                comparison_hash.update(_encoded(('left-only', a)))
                a = next(left_rows, None)
            elif a is None or b[0] < a[0]:
                events.add('extra_feature_right', key=b[0])
                comparison_hash.update(_encoded(('right-only', b)))
                b = next(right_rows, None)
            else:
                common += 1
                same = a[1] == b[1]
                equal += int(same)
                comparison_hash.update(_encoded(('pair', a, b)))
                if not same:
                    events.add('changed_feature', key=a[0], left=a[1], right=b[1])
                a, b = next(left_rows, None), next(right_rows, None)
        left_rows = _bridge_rows325(left, left_arm, expected_ranks, steps, events, 'left', bridge_stats[0])
        right_rows = _bridge_rows325(right, right_arm, expected_ranks, steps, events, 'right', bridge_stats[1])
        a, b = next(left_rows, None), next(right_rows, None)
        while a is not None or b is not None:
            if b is None or (a is not None and a[0] < b[0]):
                events.add('unpaired_bridge_left', rank=a[0][0], step=a[0][1])
                a = next(left_rows, None)
            elif a is None or b[0] < a[0]:
                events.add('unpaired_bridge_right', rank=b[0][0], step=b[0][1])
                b = next(right_rows, None)
            else:
                bridge_pairs += 1
                try:
                    result = inv.compare_inverse_bridge325(a[1], b[1])
                except (q.v.AuditError, KeyError, TypeError, ValueError, IndexError) as exc:
                    # Error text is capped; raw contexts never enter reports.
                    events.add('inverse_bridge_difference', rank=a[0][0], step=a[0][1], error=str(exc)[:512])
                    bridge_hash.update(_encoded((a[0], 'failed', type(exc).__name__, str(exc))))
                else:
                    _require(result.get('passed') is True, 'stream.bridge_verdict')
                    bridge_proved += 1
                    moved += result['moved_physical_nodes']
                    bridge_hash.update(_encoded((a[0], result)))
                a, b = next(left_rows, None), next(right_rows, None)
    def summaries(items):
        return [dict(rows=x['rows'], unique=x['unique'], sha256=x['hash'].hexdigest()) for x in items]
    failures = events.report()
    return dict(contract=CONTRACT325, passed=failures['count'] == 0,
                arms=[left_arm, right_arm], expected_ranks=expected_ranks, steps=steps,
                databases=identities, features=dict(sides=summaries(feature_stats), common=common, equal=equal,
                comparison_sha256=comparison_hash.hexdigest()), inverse_bridge=dict(sides=summaries(bridge_stats),
                expected=expected_ranks * steps, paired=bridge_pairs, proved=bridge_proved,
                moved_physical_nodes=moved, comparison_sha256=bridge_hash.hexdigest()), failures=failures)


def certify_providers325(database, arm, *, expected_ranks=256, expected_families=254,
                          expected_family_keys=None, progress=None):
    """Return compact proofs immediately after provider closure; no feature scan.

    Caller persists this result before beginning cross-arm feature comparison.
    The full original route reconciler is reused unchanged.
    """
    _require(expected_ranks == q.N, 'stream.provider_world_size', expected_ranks)
    _require(type(expected_families) is int and 0 < expected_families <= 254, 'stream.provider_family_limit')
    inv.require_inverse_source_binding325()
    with readonly_database325(database) as (connection, identity):
        identity['schema_sha256'] = _schema325(connection, ('route', 'context', 'domain'))
        _arm325(connection, ('route', 'context', 'domain'), arm)
        rank_digest = _ranks325(connection, 'domain', arm, expected_ranks)
        family_rows = connection.execute('select kind,key,count(*) from route where arm=? '
                                         'group by kind,key order by kind,key limit ?',
                                         (arm, expected_families + 1)).fetchall()
        _require(len(family_rows) == expected_families, 'stream.provider_family_census')
        census_hash = hashlib.sha256()
        families = []
        for kind, key, count in family_rows:
            _require(kind in ('routes', 'adapters', 'static') and isinstance(key, str) and
                     len(key) <= 4096 and count == expected_ranks,
                     'stream.provider_family_ranks', (kind, key, count))
            ranks = connection.execute('select rank from route where arm=? and kind=? and key=? order by rank',
                                       (arm, kind, key))
            for rank in range(expected_ranks):
                _require(next(ranks, None) == (rank,), 'stream.provider_family_rank_identity', (kind, key, rank))
            _require(next(ranks, None) is None, 'stream.provider_extra_rank', (kind, key))
            contexts = connection.execute('select rank from context where arm=? and kind=? and key=? order by rank',
                                          (arm, kind, key))
            for rank in range(expected_ranks):
                _require(next(contexts, None) == (rank,), 'stream.provider_context_rank_identity', (kind, key, rank))
            _require(next(contexts, None) is None, 'stream.provider_extra_context', (kind, key))
            families.append((kind, key))
            census_hash.update(_encoded((kind, key, count)))
        if expected_family_keys is not None:
            _require(families == sorted(tuple(x) for x in expected_family_keys), 'stream.provider_expected_families')
        extras = connection.execute("select count(*) from context c where c.arm=? and c.kind!='inverse-bridge' "
            'and not exists (select 1 from route r where r.arm=c.arm and r.kind=c.kind and r.key=c.key and r.rank=c.rank)',
            (arm,)).fetchone()[0]
        _require(extras == 0, 'stream.provider_orphan_context')
        proofs, counts = q.route_certificates(connection, arm, progress)
        _require(counts['families'] == len(proofs) == expected_families, 'stream.provider_proof_census')
        _require(set(proofs) == {repr(x) for x in families}, 'stream.provider_proof_identity')
        _require(all(isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)
                     for value in proofs.values()), 'stream.provider_proof_digest')
        return dict(contract=CONTRACT325, passed=True, arm=arm, rank_count=expected_ranks,
                    database=identity, counts=counts, proofs=proofs, proofs_sha256=_digest(proofs),
                    families_sha256=_digest(families), family_rank_census_sha256=census_hash.hexdigest(),
                    domain_rank_census_sha256=rank_digest)


def compare_provider_certificates325(left, right, *, expected_ranks=256, expected_families=254):
    """Validate compact provider certificates, then compare exact family identities."""
    _require([left.get('arm'), right.get('arm')] == ['qual-C0', 'qual-L256'],
             'stream.provider_pair_arms')
    for certificate in (left, right):
        _require(certificate.get('contract') == CONTRACT325 and certificate.get('passed') is True,
                 'stream.provider_certificate_contract')
        _require(certificate.get('rank_count') == expected_ranks and
                 certificate['counts']['families'] == len(certificate['proofs']) == expected_families,
                 'stream.provider_certificate_census')
        _require(certificate['proofs_sha256'] == _digest(certificate['proofs']),
                 'stream.provider_certificate_digest')
    events = _Events325()
    a, b = left['proofs'], right['proofs']
    for key in sorted(a.keys() | b.keys()):
        if key not in a:
            events.add('extra_provider_right', key=key)
        elif key not in b:
            events.add('missing_provider_right', key=key)
        elif a[key] != b[key]:
            events.add('changed_provider', key=key, left=a[key], right=b[key])
    failures = events.report()
    return dict(contract=CONTRACT325, passed=failures['count'] == 0,
                families=expected_families, proofs_sha256=[left['proofs_sha256'], right['proofs_sha256']],
                failures=failures)
