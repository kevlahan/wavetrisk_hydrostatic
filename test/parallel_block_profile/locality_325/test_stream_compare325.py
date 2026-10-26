"""Real SQLite/bridge/reconciler checks for bounded streaming finalization.

Constructed fixture worlds do not grant native qualification. Source binding is
controlled only for these checker tests; native callers must bind source bytes.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import tracemalloc
import unittest

import inverse_identity325 as inv
import stream_compare325 as s
import verify325 as q
from test_inverse_bridge_adversary325 import schedule
from test_provider_paths325 import inverses


def digest(value):
    return hashlib.sha256(str(value).encode()).hexdigest()


def fixture(path, arm, *, ranks=2, steps=2, features=8, duplicates=False):
    inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})
    scopes, plan = schedule((0, 1, 2) if arm == 'qual-C0' else (0, 2, 1))
    context = inv.inverse_bridge_context325(scopes, plan)
    c = q._create_cache(path)
    try:
        if duplicates:
            c.execute('drop table feature')
            c.execute('create table feature(arm text,key text,digest text)')
            c.execute('drop table context')
            c.execute('create table context(arm text,rank int,kind text,key text,payload text)')
        c.executemany('insert into domain values(?,?,?)', ((arm, rank, '[]') for rank in range(ranks)))
        c.executemany('insert into feature values(?,?,?)',
                      ((arm, 'key-' + str(i).zfill(8), digest(i)) for i in range(features)))
        for rank in range(ranks):
            for step in range(1, steps + 1):
                context['step'] = context['epoch'] = step
                c.execute('insert into context values(?,?,?,?,?)',
                          (arm, rank, 'inverse-bridge', str(step), json.dumps(context)))
        c.commit()
    finally:
        c.close()
    return path


def change(path, statement, params=()):
    c = sqlite3.connect(path)
    try:
        c.execute(statement, params)
        c.commit()
    finally:
        c.close()


class StreamingComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.left = fixture(self.root / 'left.sqlite', 'qual-C0')
        self.right = fixture(self.root / 'right.sqlite', 'qual-L256')

    def compare(self, **kwargs):
        return s.compare_databases325(self.left, self.right, expected_ranks=2, steps=2, **kwargs)

    def test_exact_streaming_pair_and_real_physical_bridge(self):
        before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in (self.left, self.right)]
        result = self.compare()
        self.assertTrue(result['passed'])
        self.assertEqual(result['features']['equal'], 8)
        self.assertEqual(result['inverse_bridge']['proved'], 4)
        self.assertEqual(result['inverse_bridge']['moved_physical_nodes'], 8)
        self.assertEqual(before, [hashlib.sha256(path.read_bytes()).hexdigest() for path in (self.left, self.right)])
        self.assertNotIn('features', result['features'])
        self.assertLess(len(json.dumps(result)), 6000)

    def test_changed_missing_and_extra_feature_each_fail(self):
        change(self.right, "delete from feature where key='key-00000000'")
        change(self.right, "update feature set digest=? where key='key-00000001'", (digest('wrong'),))
        change(self.right, 'insert into feature values(?,?,?)', ('qual-L256', 'extra', digest('extra')))
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['counts'],
                         dict(missing_feature_right=1, extra_feature_right=1, changed_feature=1))
        self.assertEqual(result['features']['sides'][0]['rows'], 8)
        self.assertEqual(result['features']['common'], 7)

    def test_duplicate_features_even_if_equal_digest_fail(self):
        self.left.unlink()
        fixture(self.left, 'qual-C0', duplicates=True)
        change(self.left, "insert into feature select arm,key,digest from feature where key='key-00000000'")
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['counts']['duplicate_feature'], 1)

    def test_empty_feature_pair_cannot_pass(self):
        for path in (self.left, self.right):
            change(path, 'delete from feature')
        result = self.compare()
        self.assertEqual(result['failures']['counts']['empty_features'], 2)
        self.assertFalse(result['passed'])

    def test_complete_mismatch_census_has_bounded_examples_and_repeatable_hash(self):
        change(self.right, 'update feature set digest=?', (digest('changed'),))
        result = self.compare()
        self.assertEqual(result['failures']['count'], 8)
        self.assertEqual(result['failures']['sha256'], self.compare()['failures']['sha256'])
        for i in range(100):
            change(self.right, 'insert into feature values(?,?,?)', ('qual-L256', 'extra-' + str(i), digest(i)))
        result = self.compare()
        self.assertEqual(result['failures']['count'], 108)
        self.assertEqual(len(result['failures']['examples']), 32)
        self.assertLess(len(json.dumps(result)), 16000)

    def test_long_corrupt_keys_cannot_expand_example_payload_without_limit(self):
        change(self.right, 'insert into feature values(?,?,?)', ('qual-L256', 'x' * 100000, digest('x')))
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertLess(len(json.dumps(result)), 7000)

    def test_missing_same_bridge_in_both_worlds_cannot_pass(self):
        for path in (self.left, self.right):
            change(path, "delete from context where rank=1 and key='2'")
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['counts']['missing_bridge'], 2)
        self.assertEqual(result['inverse_bridge']['paired'], 3)

    def test_missing_single_bridge_cannot_pass(self):
        change(self.right, "delete from context where rank=0 and key='2'")
        result = self.compare()
        self.assertEqual(result['failures']['counts'], dict(missing_bridge=1, unpaired_bridge_left=1))

    def test_duplicate_bridge_even_identical_payload_fails(self):
        self.right.unlink()
        fixture(self.right, 'qual-L256', duplicates=True)
        change(self.right, "insert into context select * from context where rank=0 and key='1'")
        result = self.compare()
        self.assertEqual(result['failures']['counts']['duplicate_bridge'], 1)
        self.assertFalse(result['passed'])

    def test_extra_rank_and_noncanonical_step_are_structural_failures(self):
        change(self.right, "update context set rank=2 where rank=1")
        with self.assertRaisesRegex(q.v.AuditError, 'bridge_key'):
            self.compare()
        change(self.right, 'update context set rank=1 where rank=2')
        change(self.right, "update context set key='01' where key='1'")
        with self.assertRaisesRegex(q.v.AuditError, 'bridge_key'):
            self.compare()

    def test_bridge_row_step_must_match_context(self):
        c = sqlite3.connect(self.right)
        value = json.loads(c.execute("select payload from context where rank=0 and key='1'").fetchone()[0])
        c.close()
        value['step'] = 2
        change(self.right, "update context set payload=? where rank=0 and key='1'", (json.dumps(value),))
        with self.assertRaisesRegex(q.v.AuditError, 'bridge_step_binding'):
            self.compare()

    def test_bridge_physical_byte_and_axis_corruption_are_rejected(self):
        c = sqlite3.connect(self.right)
        value = json.loads(c.execute("select payload from context where rank=0 and key='1'").fetchone()[0])
        c.close()
        value['publication']['(1, 1, 1)'][-1][-1] = 'wrong physical bytes'
        change(self.right, "update context set payload=? where rank=0 and key='1'", (json.dumps(value),))
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['counts']['inverse_bridge_difference'], 1)
        self.assertIn('physical_publication_difference', result['failures']['examples'][0]['error'])

    def test_changed_raw_bridge_values_and_auxiliary_operation_axis_fail(self):
        c = sqlite3.connect(self.right)
        value = json.loads(c.execute("select payload from context where rank=0 and key='1'").fetchone()[0])
        c.close()
        value['phases']['2']['inverse.value']['digest'] = digest('changed')
        change(self.right, "update context set payload=? where rank=0 and key='1'", (json.dumps(value),))
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertIn('raw_value_difference', result['failures']['examples'][0]['error'])
        value['bindings']['inverse.scaffold'] = [[101, 1]]
        change(self.right, "update context set payload=? where rank=0 and key='1'", (json.dumps(value),))
        result = self.compare()
        self.assertFalse(result['passed'])
        self.assertIn('auxiliary', result['failures']['examples'][0]['error'])

    def test_all_bridge_errors_counted_with_only32_examples(self):
        for path in (self.left, self.right):
            path.unlink()
        fixture(self.left, 'qual-C0', ranks=20, steps=2)
        fixture(self.right, 'qual-L256', ranks=20, steps=2)
        c = sqlite3.connect(self.right)
        for rank, key, payload in list(c.execute('select rank,key,payload from context')):
            value = json.loads(payload)
            value['phases']['3']['inverse.value']['digest'] = digest('wrong')
            c.execute('update context set payload=? where rank=? and key=?', (json.dumps(value), rank, key))
        c.commit()
        c.close()
        result = s.compare_databases325(self.left, self.right, expected_ranks=20, steps=2)
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['counts'], dict(inverse_bridge_difference=40))
        self.assertEqual(len(result['failures']['examples']), 32)
        self.assertEqual(result['inverse_bridge']['paired'], 40)
        self.assertEqual(result['inverse_bridge']['proved'], 0)
        self.assertLess(len(json.dumps(result)), 16000)

    def test_bridge_missing_source_binding_cannot_pass(self):
        with self.assertRaises(q.v.AuditError):
            inv.bind_inverse_bridge_sources325({'candidate': {}})
        with self.assertRaisesRegex(q.v.AuditError, 'unbound_source'):
            self.compare()

    def test_wrong_arm_and_incomplete_domain_cannot_pass(self):
        change(self.right, "update domain set arm='qual-C0'")
        with self.assertRaisesRegex(q.v.AuditError, 'unexpected_arm'):
            self.compare()
        change(self.right, "update domain set arm='qual-L256'")
        change(self.right, 'delete from domain where rank=1')
        with self.assertRaisesRegex(q.v.AuditError, 'rank_census'):
            self.compare()

    def test_schema_views_and_changed_columns_rejected(self):
        change(self.right, 'alter table feature rename column digest to hash')
        with self.assertRaisesRegex(q.v.AuditError, 'schema_columns'):
            self.compare()

    def test_real_sorted_cursor_does_not_allocate_full_feature_dictionary(self):
        for path, arm in ((self.left, 'qual-C0'), (self.right, 'qual-L256')):
            change(path, 'delete from feature')
            c = sqlite3.connect(path)
            c.executemany('insert into feature values(?,?,?)',
                          ((arm, 'key-' + str(i).zfill(8), digest(i)) for i in range(25000)))
            c.commit()
            c.close()
        tracemalloc.start()
        result = self.compare()
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.assertTrue(result['passed'])
        self.assertEqual(result['features']['equal'], 25000)
        self.assertLess(peak, 3 * 1024 * 1024)
        self.assertLess(len(json.dumps(result)), 6000)

    def test_original_and_streamed_comparators_agree_on_full_census(self):
        # This uses the real original comparator with all256x4bridge rows.
        for path in (self.left, self.right):
            path.unlink()
        fixture(self.left, 'qual-C0', ranks=256, steps=4)
        fixture(self.right, 'qual-L256', ranks=256, steps=4)
        c = q._create_cache(':memory:')
        try:
            for path in (self.left, self.right):
                source = sqlite3.connect(path)
                for table in ('domain', 'feature', 'context'):
                    rows = source.execute('select * from ' + table)
                    placeholders = {'domain': '?,?,?', 'feature': '?,?,?', 'context': '?,?,?,?,?'}[table]
                    c.executemany('insert into ' + table + ' values(' + placeholders + ')', rows)
                source.close()
            self.assertTrue(q.compare_features(c, 'qual-C0', 'qual-L256')['passed'])
            self.assertTrue(s.compare_databases325(self.left, self.right)['passed'])
            c.execute("update feature set digest=? where arm='qual-L256'", (digest('wrong'),))
            change(self.right, 'update feature set digest=?', (digest('wrong'),))
            with self.assertRaisesRegex(q.v.AuditError, 'ordered_program_or_value_difference'):
                q.compare_features(c, 'qual-C0', 'qual-L256')
            self.assertFalse(s.compare_databases325(self.left, self.right)['passed'])
        finally:
            c.close()


class ReadonlyDatabaseTests(unittest.TestCase):
    def test_readonly_connection_rejects_write(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder).resolve() / 'db.sqlite'
            c = sqlite3.connect(path)
            c.execute('create table foo(value)')
            c.close()
            before = path.read_bytes()
            with s.readonly_database325(path) as (c, _):
                with self.assertRaises(sqlite3.OperationalError):
                    c.execute('insert into foo values(1)')
            self.assertEqual(before, path.read_bytes())

    def test_uri_special_characters_are_literal_path(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder).resolve() / 'a?mode=rw#b.sqlite'
            c = sqlite3.connect(path)
            c.execute('create table foo(value)')
            c.close()
            with s.readonly_database325(path) as (c, identity):
                self.assertEqual(identity['path'], str(path))
                self.assertEqual(c.execute('pragma query_only').fetchone()[0], 1)

    def test_symlink_file_parent_and_sidecars_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder).resolve()
            actual = root / 'actual'
            actual.mkdir()
            db = actual / 'db.sqlite'
            sqlite3.connect(db).close()
            link = root / 'db.sqlite'
            link.symlink_to(db)
            directory = root / 'directory'
            directory.symlink_to(actual)
            for path in (link, directory / 'db.sqlite'):
                with self.subTest(path=path), self.assertRaisesRegex(q.v.AuditError, 'symlink'):
                    with s.readonly_database325(path):
                        pass
            for suffix in ('-wal', '-shm', '-journal'):
                sidecar = Path(str(db) + suffix)
                sidecar.write_bytes(b'')
                with self.subTest(suffix=suffix), self.assertRaisesRegex(q.v.AuditError, 'sidecar'):
                    with s.readonly_database325(db):
                        pass
                sidecar.unlink()

    def test_wal_header_rejected_even_with_no_transaction_sidecars(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder).resolve() / 'db.sqlite'
            c = sqlite3.connect(path)
            c.execute('pragma journal_mode=wal')
            c.execute('create table foo(value)')
            c.commit()
            c.close()
            with self.assertRaisesRegex(q.v.AuditError, 'wal_mode'):
                with s.readonly_database325(path):
                    pass


class ProviderCertificatesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name).resolve()
        cls.database = cls.root / 'providers.sqlite'
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})
        plans, contexts, _ = inverses('inverse.alias', unallocated=True)
        c = q._create_cache(cls.database)
        for rank in range(256):
            c.execute('insert into domain values(?,?,?)', ('qual-C0', rank, json.dumps([rank])))
        for family in range(254):
            key = 'inverse.alias/' + str(family).zfill(3)
            for rank in range(256):
                c.execute('insert into route values(?,?,?,?,?)',
                          ('qual-C0', rank, 'static', key, json.dumps(plans[rank])))
                c.execute('insert into context values(?,?,?,?,?)',
                          ('qual-C0', rank, 'static', key, json.dumps(contexts[rank])))
        c.commit()
        c.close()
        cls.certificate = s.certify_providers325(cls.database, 'qual-C0')

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def setUp(self):
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})

    def test_real_reconcilers_complete254_without_feature_table_read(self):
        self.assertTrue(self.certificate['passed'])
        self.assertEqual(len(self.certificate['proofs']), 254)
        self.assertEqual(self.certificate['counts']['families'], 254)
        self.assertNotIn('features', self.certificate)
        self.assertNotIn('bridges', self.certificate)
        self.assertLess(len(json.dumps(self.certificate)), 40000)
        other = copy.deepcopy(self.certificate)
        other['arm'] = 'qual-L256'
        self.assertTrue(s.compare_provider_certificates325(self.certificate, other)['passed'])

    def test_wrong_family_and_changed_proof_cannot_pass(self):
        other = copy.deepcopy(self.certificate)
        other['arm'] = 'qual-L256'
        key = next(iter(other['proofs']))
        other['proofs']['wrong-family'] = other['proofs'].pop(key)
        other['proofs_sha256'] = s._digest(other['proofs'])
        result = s.compare_provider_certificates325(self.certificate, other)
        self.assertFalse(result['passed'])
        self.assertEqual(result['failures']['count'], 2)
        other = copy.deepcopy(self.certificate)
        other['arm'] = 'qual-L256'
        other['proofs'][key] = digest('wrong')
        other['proofs_sha256'] = s._digest(other['proofs'])
        self.assertEqual(s.compare_provider_certificates325(self.certificate, other)['failures']['counts'],
                         dict(changed_provider=1))

    def test_wrong_provider_pair_identity_is_rejected(self):
        with self.assertRaisesRegex(q.v.AuditError, 'pair_arms'):
            s.compare_provider_certificates325(self.certificate, self.certificate)

    def test_tampered_proof_hash_or_missing_family_fails_closed(self):
        other = copy.deepcopy(self.certificate)
        other['arm'] = 'qual-L256'
        other['proofs_sha256'] = 'f' * 64
        with self.assertRaisesRegex(q.v.AuditError, 'certificate_digest'):
            s.compare_provider_certificates325(self.certificate, other)
        del other['proofs'][next(iter(other['proofs']))]
        with self.assertRaisesRegex(q.v.AuditError, 'certificate_census'):
            s.compare_provider_certificates325(self.certificate, other)

    def test_exact_expected_family_identity_is_enforced(self):
        keys = [('static', 'wrong/' + str(i)) for i in range(254)]
        with self.assertRaisesRegex(q.v.AuditError, 'expected_families'):
            s.certify_providers325(self.database, 'qual-C0', expected_family_keys=keys)

    def test_missing_route_rank_even_with254families_is_rejected(self):
        other = self.root / 'bad-rank.sqlite'
        other.write_bytes(self.database.read_bytes())
        self.addCleanup(lambda: other.unlink(missing_ok=True))
        change(other, "delete from route where rank=255 and key='inverse.alias/000'")
        with self.assertRaisesRegex(q.v.AuditError, 'family_ranks'):
            s.certify_providers325(other, 'qual-C0')

    def test_missing_context_and_orphan_context_are_rejected(self):
        other = self.root / 'bad-context.sqlite'
        other.write_bytes(self.database.read_bytes())
        self.addCleanup(lambda: other.unlink(missing_ok=True))
        change(other, "delete from context where rank=255 and key='inverse.alias/000'")
        with self.assertRaisesRegex(q.v.AuditError, 'context_rank_identity'):
            s.certify_providers325(other, 'qual-C0')
        other.write_bytes(self.database.read_bytes())
        change(other, "insert into context values('qual-C0',0,'static','orphan','{}')")
        with self.assertRaisesRegex(q.v.AuditError, 'orphan_context'):
            s.certify_providers325(other, 'qual-C0')


if __name__ == '__main__':
    unittest.main()
