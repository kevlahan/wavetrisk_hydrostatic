"""Deterministic, bounded per-rank cache construction for saved Stage325 traces.

This module changes execution scheduling only. Each worker calls the unchanged
cache_rank325 gate and hashes the actual compressed bytes through the recovery
reader. SQLite files are never shared by writers. A merged file is usable only
after every shard's exact rows, primary keys and table digests reconcile.
"""
import hashlib
import importlib
import json
import multiprocessing
from pathlib import Path
import sqlite3

import verify325 as q

TABLES = ('route', 'domain', 'context', 'feature')
KEYS = {'route': ('arm', 'rank', 'kind', 'key'),
        'domain': ('arm', 'rank'),
        'context': ('arm', 'rank', 'kind', 'key'),
        'feature': ('arm', 'key')}
ARMS = ('qual-C0', 'qual-L256')


def require(ok, message):
    q.require(ok, 'parallel325.' + message)


def scratch_size325(root, limit=q.SCRATCH_LIMIT):
    size = sum(p.stat().st_size for p in Path(root).rglob('*') if p.is_file())
    require(size <= limit, 'aggregate_scratch_bound')
    return size


def _regular_inside(path, root):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'regular_shard_required')
    require(path.resolve().is_relative_to(Path(root).resolve()), 'shard_outside_scratch')
    return path.resolve()


def _file_hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def table_manifest325(connection, table, *, source=None):
    """Exact ordered row digest, optionally for rows selected by a joined shard."""
    require(table in TABLES, 'unknown_table')
    keys = KEYS[table]
    if source is None:
        statement = 'select * from ' + table + ' order by ' + ','.join(keys)
    else:
        require(source == 'shard', 'unknown_join')
        match = ' and '.join(f'd.{k}=s.{k}' for k in keys)
        order = ','.join('s.' + k for k in keys)
        statement = (f'select d.* from main.{table} d join shard.{table} s on '
                     + match + ' order by ' + order)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(statement):
        digest.update((json.dumps(row, separators=(',', ':'), ensure_ascii=True) + '\n').encode('ascii'))
        count += 1
    return {'rows': count, 'sha256': digest.hexdigest()}


def shard_certificate325(path):
    path = Path(path)
    connection = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)
    try:
        tables = {t: table_manifest325(connection, t) for t in TABLES}
    finally:
        connection.close()
    return {'path': str(path), 'sha256': _file_hash(path), 'tables': tables}


def merge_shards325(shards, destination, root, *, limit=q.SCRATCH_LIMIT):
    """Merge exact disjoint rows; remove a shard only after verified commit.

    Each input is a shard_certificate325 result. On any error the destination
    is partial evidence, never a qualified cache. Fresh scratch files alone use
    journal_mode=OFF, avoiding temporary rollback copies that defeat the cap.
    """
    root = Path(root).resolve(strict=True)
    destination = Path(destination)
    require(destination.parent.resolve() == root and not destination.exists(),
            'fresh_destination_required')
    require(bool(shards), 'no_shards')
    paths = [_regular_inside(s['path'], root) for s in shards]
    require(len(set(paths)) == len(paths), 'duplicate_shard')
    # URI interpretation is a per-connection option on bbserv's SQLite build.
    # Create the checked schema, then explicitly reopen with URI support before
    # attaching read-only file: shard names. Never rely on a compile-time default.
    connection = q._create_cache(destination)
    connection.close()
    connection = sqlite3.connect(destination.resolve().as_uri(), uri=True)
    connection.execute('pragma journal_mode=OFF')
    counts = dict.fromkeys(TABLES, 0)
    proofs = []
    peak = scratch_size325(root, limit)
    try:
        for certificate, path in zip(shards, paths):
            require(_file_hash(path) == certificate['sha256'], 'shard_file_changed')
            # Every page allocated to main must fit beside all remaining shards.
            total = scratch_size325(root, limit)
            page = connection.execute('pragma page_size').fetchone()[0]
            current = destination.stat().st_size
            max_pages = (limit - (total - current)) // page
            connection.execute('pragma max_page_count=' + str(max_pages))
            connection.execute('attach database ? as shard', (path.as_uri() + '?mode=ro',))
            try:
                source = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)
                try:
                    actual = {t: table_manifest325(source, t) for t in TABLES}
                finally:
                    source.close()
                require(actual == certificate['tables'], 'shard_table_changed')
                for table in TABLES:
                    # Plain INSERT makes every duplicated primary key fatal.
                    connection.execute(f'insert into main.{table} select * from shard.{table}')
                connection.commit()
                joined = {t: table_manifest325(connection, t, source='shard') for t in TABLES}
                require(joined == actual, 'merged_rows_changed')
                for table in TABLES:
                    counts[table] += actual[table]['rows']
                    actual_count = connection.execute('select count(*) from ' + table).fetchone()[0]
                    require(actual_count == counts[table], 'merged_row_census')
                peak = max(peak, scratch_size325(root, limit))
                proofs.append({'source_sha256': certificate['sha256'], 'tables': joined})
            finally:
                # Release an aborted INSERT transaction before DETACH so the
                # original constraint/corruption error is never masked.
                connection.rollback()
                connection.execute('detach database shard')
            path.unlink()
            peak = max(peak, scratch_size325(root, limit))
    finally:
        connection.close()
    return {'database': str(destination.resolve()), 'row_counts': counts,
            'verified_shards': proofs, 'sampled_peak_scratch_bytes': peak,
            'database_sha256': _file_hash(destination)}


def bind_worker_sources325(recovery, original):
    """Independently bind every spawned worker to the original compiled code."""
    binding_path = Path(recovery.__file__).resolve().parent / 'binding325.json'
    require(binding_path.is_file() and not binding_path.is_symlink(),
            'worker_binding_required')
    binding = json.loads(binding_path.read_text())
    original = Path(original).resolve(strict=True)
    require(original == Path(binding['original_results']).resolve(strict=True),
            'worker_original_results_binding')
    manifest_path = Path(binding['original_kit']) / 'source-files.json'
    require(manifest_path.is_file() and not manifest_path.is_symlink(),
            'worker_source_manifest_required')
    require(_file_hash(manifest_path) == binding['source_manifest_sha256'],
            'worker_source_manifest_changed')
    manifest = json.loads(manifest_path.read_text())
    from inverse_identity325 import bind_inverse_bridge_sources325
    proof = bind_inverse_bridge_sources325(
        manifest, source_root=original / 'sources' / 'candidate')
    require(proof.get('passed') is True, 'worker_inverse_source_binding')
    return {**proof, 'source_manifest_sha256': binding['source_manifest_sha256']}


def _worker325(job):
    arm, index, tasks, database, root, original, anchors, limit, quota, reader_module = job
    observed = {}
    recovery = importlib.import_module(reader_module)
    source_binding = bind_worker_sources325(recovery, original)
    original_reader = q.v.open_stream
    q.v.open_stream = lambda path: recovery.verified_stream(path, original, anchors, observed)
    connection = q._create_cache(database)
    connection.execute('pragma journal_mode=OFF')
    page = connection.execute('pragma page_size').fetchone()[0]
    connection.execute('pragma max_page_count=' + str(quota // page))
    summaries = []
    peak = scratch_size325(root, limit)
    try:
        for rank, path in tasks:
            summaries.append(q.cache_rank325(connection, path, arm, rank))
            peak = max(peak, scratch_size325(root, limit))
    finally:
        connection.close()
        q.v.open_stream = original_reader
    require(len(observed) == len(tasks), 'worker_actual_reads_incomplete')
    certificate = shard_certificate325(Path(database).resolve())
    return {'arm': arm, 'index': index, 'certificate': certificate,
            'streams': summaries, 'capture_reads': observed,
            'source_binding': source_binding,
            'sampled_peak_scratch_bytes': peak}


def tasks_from_results325(results):
    results = Path(results)
    return [{'arm': arm, 'rank': rank, 'path': str(path.resolve())}
            for arm in ARMS for rank, path in enumerate(q.v._stream_set(results / arm, q.N))]


def cache_shards325(tasks, root, workers=16, *, original_root, anchors,
                    limit=q.SCRATCH_LIMIT, reader_module='recover325', progress=None):
    """Build two complete optimized caches using rank-strided private writers.

    Caller owns the allocation/alarm/process-group and aggregate-memory guard.
    No source, solver or predecessor artifact is written by this function.
    """
    root = Path(root)
    require(root.is_dir() and not any(root.iterdir()), 'fresh_scratch_required')
    require(isinstance(workers, int) and 2 <= workers <= 16 and workers % 2 == 0,
            'worker_count')
    require(len(tasks) == 2 * q.N, 'task_count')
    by_arm = {arm: {} for arm in ARMS}
    for task in tasks:
        arm, rank, path = task['arm'], task['rank'], Path(task['path'])
        require(arm in ARMS and isinstance(rank, int) and 0 <= rank < q.N,
                'task_identity')
        require(rank not in by_arm[arm], 'duplicate_rank')
        require(path.is_file() and not path.is_symlink(), 'regular_capture_required')
        require(path.resolve().is_relative_to(Path(original_root).resolve()), 'capture_outside_original')
        by_arm[arm][rank] = str(path.resolve())
    require(all(sorted(rows) == list(range(q.N)) for rows in by_arm.values()), 'missing_rank')
    jobs = []
    n_shards = workers // 2
    quota = limit // (workers + 1)  # Keep one shard's space free for the merge.
    for arm in ARMS:
        for index in range(n_shards):
            selected = [(rank, by_arm[arm][rank]) for rank in range(index, q.N, n_shards)]
            database = root / f'{arm}-shard-{index}.sqlite'
            jobs.append((arm, index, selected, str(database), str(root), str(original_root),
                         anchors, limit, quota, reader_module))
    finished = {arm: [] for arm in ARMS}
    peak = 0
    # Spawn guarantees every worker opens its own SQLite connection and reader.
    context = multiprocessing.get_context('spawn')
    # multiprocessing.Pool's context exit terminates outstanding workers on
    # an exception/alarm. Executor.__exit__ would wait for their full traces.
    with context.Pool(processes=workers) as pool:
        for result in pool.imap_unordered(_worker325, jobs):
            finished[result['arm']].append(result)
            peak = max(peak, result['sampled_peak_scratch_bytes'], scratch_size325(root, limit))
            if progress:
                progress('shard', result['arm'], result['index'], len(result['streams']))
    report = {'databases': {}, 'streams': {}, 'capture_reads': {}, 'merges': {},
              'source_bindings': {},
              'workers': workers, 'sampled_peak_scratch_bytes': peak}
    for arm in ARMS:
        shards = sorted(finished[arm], key=lambda row: row['index'])
        require([row['index'] for row in shards] == list(range(n_shards)), 'missing_shard')
        streams = sorted((s for row in shards for s in row['streams']), key=lambda row: row['rank'])
        require([s['rank'] for s in streams] == list(range(q.N)), 'stream_rank_census')
        for row in shards:
            for name, check in row['capture_reads'].items():
                require(name not in report['capture_reads'], 'duplicated_actual_read')
                report['capture_reads'][name] = check
        merged = merge_shards325([row['certificate'] for row in shards], root / (arm + '.sqlite'), root, limit=limit)
        report['merges'][arm] = merged
        report['databases'][arm] = merged['database']
        report['streams'][arm] = streams
        report['source_bindings'][arm] = [row['source_binding'] for row in shards]
        report['sampled_peak_scratch_bytes'] = max(report['sampled_peak_scratch_bytes'], merged['sampled_peak_scratch_bytes'])
    require(len(report['capture_reads']) == 2 * q.N, 'actual_read_census')
    return report
