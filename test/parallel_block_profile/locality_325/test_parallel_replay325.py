"""Exact SQLite merge tests for parallel replay scheduling; no native verdict."""
import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace
import unittest

import parallel_replay325 as p
import verify325 as q


def fixture(path, rank):
    connection = q._create_cache(path)
    connection.execute('insert into route values(?,?,?,?,?)', ('qual-C0', rank, 'routes', 'same-key', json.dumps({'rank': rank, 'wire': [3, 2, 1]})))
    connection.execute('insert into domain values(?,?,?)', ('qual-C0', rank, json.dumps([rank])))
    connection.execute('insert into context values(?,?,?,?,?)', ('qual-C0', rank, 'routes', 'same-key', json.dumps({'rows': [rank]})))
    connection.execute('insert into feature values(?,?,?)', ('qual-C0', 'physical-' + str(rank), 'value-' + str(rank)))
    connection.commit()
    connection.close()
    return p.shard_certificate325(path)


class ParallelMergeTests(unittest.TestCase):
    def test_exact_serial_and_merged_tables(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            shards = [fixture(root / f'shard-{rank}.sqlite', rank) for rank in (2, 0, 1)]
            result = p.merge_shards325(shards, root / 'merged.sqlite', root)
            self.assertEqual(result['row_counts'], dict.fromkeys(p.TABLES, 3))
            self.assertEqual(len(result['verified_shards']), 3)
            self.assertTrue(all(not Path(s['path']).exists() for s in shards))
            merged = sqlite3.connect(root / 'merged.sqlite')
            reference = q._create_cache(':memory:')
            try:
                for rank in (0, 1, 2):
                    source = sqlite3.connect(fixture(root / f'reference-{rank}.sqlite', rank)['path'])
                    for table in p.TABLES:
                        rows = source.execute('select * from ' + table).fetchall()
                        reference.executemany('insert into ' + table + ' values(' + ','.join('?' for _ in rows[0]) + ')', rows)
                    source.close()
                for table in p.TABLES:
                    self.assertEqual(p.table_manifest325(merged, table), p.table_manifest325(reference, table))
            finally:
                merged.close()
                reference.close()

    def test_changed_shard_rejected_and_retained(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            certificate = fixture(root / 'one.sqlite', 0)
            connection = sqlite3.connect(certificate['path'])
            connection.execute("update feature set digest='corrupted'")
            connection.commit()
            connection.close()
            with self.assertRaisesRegex(q.v.AuditError, 'shard_file_changed'):
                p.merge_shards325([certificate], root / 'merged.sqlite', root)
            self.assertTrue(Path(certificate['path']).exists())

    def test_changed_table_certificate_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            certificate = fixture(root / 'one.sqlite', 0)
            certificate['tables']['feature']['sha256'] = '0' * 64
            with self.assertRaisesRegex(q.v.AuditError, 'shard_table_changed'):
                p.merge_shards325([certificate], root / 'merged.sqlite', root)
            self.assertTrue(Path(certificate['path']).exists())

    def test_duplicate_primary_key_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            shards = [fixture(root / f'shard-{i}.sqlite', 0) for i in (0, 1)]
            with self.assertRaises(sqlite3.IntegrityError):
                p.merge_shards325(shards, root / 'merged.sqlite', root)
            self.assertTrue(Path(shards[1]['path']).exists())

    def test_repeated_shard_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            certificate = fixture(root / 'one.sqlite', 0)
            with self.assertRaisesRegex(q.v.AuditError, 'duplicate_shard'):
                p.merge_shards325([certificate, certificate], root / 'merged.sqlite', root)

    def test_shard_outside_fresh_scratch_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inside = root / 'scratch'
            inside.mkdir()
            certificate = fixture(root / 'outside.sqlite', 0)
            with self.assertRaisesRegex(q.v.AuditError, 'shard_outside_scratch'):
                p.merge_shards325([certificate], inside / 'merged.sqlite', inside)

    def test_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            certificate = fixture(root / 'one.sqlite', 0)
            link = root / 'link.sqlite'
            link.symlink_to(root / 'one.sqlite')
            certificate['path'] = str(link)
            with self.assertRaisesRegex(q.v.AuditError, 'regular_shard_required'):
                p.merge_shards325([certificate], root / 'merged.sqlite', root)

    def test_aggregate_limit_counts_all_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'first').write_bytes(b'123')
            (root / 'second').write_bytes(b'456')
            with self.assertRaisesRegex(q.v.AuditError, 'aggregate_scratch_bound'):
                p.scratch_size325(root, limit=5)

    def test_task_count_fails_before_process_creation(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(q.v.AuditError, 'task_count'):
                p.cache_shards325([], folder, original_root=folder, anchors={})

    def test_worker_count_bounded(self):
        with tempfile.TemporaryDirectory() as folder:
            for count in (0, 1, 3, 18):
                with self.assertRaisesRegex(q.v.AuditError, 'worker_count'):
                    p.cache_shards325([], folder, count, original_root=folder, anchors={})

    def test_duplicate_rank_fails_before_process_creation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            captures = root / 'captures'
            scratch = root / 'scratch'
            captures.mkdir()
            scratch.mkdir()
            capture = captures / 'trace.gz'
            capture.write_bytes(b'only path validation is reached')
            task = {'arm': 'qual-C0', 'rank': 0, 'path': str(capture)}
            with self.assertRaisesRegex(q.v.AuditError, 'duplicate_rank'):
                p.cache_shards325([task] * 512, scratch, original_root=captures, anchors={})

    def test_full_provider_proof_matches_serial_world(self):
        from test_provider_paths325 import routes
        plans, contexts = routes()
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            paths = [root / f'part-{i}.sqlite' for i in (0, 1)]
            shards = [q._create_cache(path) for path in paths]
            serial = q._create_cache(':memory:')
            try:
                for rank, (plan, context) in enumerate(zip(plans, contexts)):
                    for connection in (serial, shards[rank % 2]):
                        connection.execute('insert into route values(?,?,?,?,?)', ('qual-C0', rank, 'routes', 'world', json.dumps(plan)))
                        connection.execute('insert into context values(?,?,?,?,?)', ('qual-C0', rank, 'routes', 'world', json.dumps(context)))
                        connection.execute('insert into domain values(?,?,?)', ('qual-C0', rank, json.dumps([])))
                for connection in shards:
                    connection.commit()
                    connection.close()
                certificates = [p.shard_certificate325(path) for path in paths]
                merged = p.merge_shards325(certificates, root / 'merged.sqlite', root)
                connection = sqlite3.connect(merged['database'])
                try:
                    self.assertEqual(q.route_certificates(serial, 'qual-C0'), q.route_certificates(connection, 'qual-C0'))
                finally:
                    connection.close()
            finally:
                serial.close()

    def test_worker_requires_sealed_binding(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            recovery = SimpleNamespace(__file__=str(root / 'recover325.py'))
            with self.assertRaisesRegex(q.v.AuditError, 'worker_binding_required'):
                p.bind_worker_sources325(recovery, root)

    def test_worker_rejects_rebound_original_root(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            original = root / 'original'
            original.mkdir()
            (root / 'binding325.json').write_text(json.dumps({'original_results': str(original)}))
            recovery = SimpleNamespace(__file__=str(root / 'recover325.py'))
            with self.assertRaisesRegex(q.v.AuditError, 'worker_original_results_binding'):
                p.bind_worker_sources325(recovery, root)

    def test_worker_rejects_changed_manifest(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            oldkit = root / 'oldkit'
            oldkit.mkdir()
            (oldkit / 'source-files.json').write_text('{}')
            (root / 'binding325.json').write_text(json.dumps({
                'original_results': str(root), 'original_kit': str(oldkit),
                'source_manifest_sha256': '0' * 64}))
            recovery = SimpleNamespace(__file__=str(root / 'recover325.py'))
            with self.assertRaisesRegex(q.v.AuditError, 'worker_source_manifest_changed'):
                p.bind_worker_sources325(recovery, root)

    def test_worker_binds_actual_expected_source_bytes(self):
        from inverse_identity325 import EXPECTED_SOURCE_HASHES325
        project = Path(__file__).resolve().parent / 'fixture/source-binding'
        if not project.is_dir():project = Path(__file__).resolve().parents[3]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            oldkit = root / 'oldkit'
            oldkit.mkdir()
            for relative in EXPECTED_SOURCE_HASHES325:
                target = root / 'sources' / 'candidate' / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes((project / relative).read_bytes())
            manifest = oldkit / 'source-files.json'
            manifest.write_text(json.dumps({'candidate': EXPECTED_SOURCE_HASHES325}))
            (root / 'binding325.json').write_text(json.dumps({
                'original_results': str(root), 'original_kit': str(oldkit),
                'source_manifest_sha256': p._file_hash(manifest)}))
            recovery = SimpleNamespace(__file__=str(root / 'recover325.py'))
            result = p.bind_worker_sources325(recovery, root)
            self.assertTrue(result['passed'])
            self.assertEqual(result['sources'], EXPECTED_SOURCE_HASHES325)
            (root / 'sources' / 'candidate' / 'src/domain.f90').write_text('corrupt')
            with self.assertRaisesRegex(q.v.AuditError, 'source_bytes'):
                p.bind_worker_sources325(recovery, root)


if __name__ == '__main__':
    unittest.main()
