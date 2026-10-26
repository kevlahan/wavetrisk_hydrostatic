"""Offline tests of native adapter guards; no original solver/capture is run."""
import json
import copy
import multiprocessing
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest

import resume_native325 as native


def _identity325(value):
    return value


class NativeAdapter325Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.module_root = Path(__file__).resolve().parent
        self.frozen = self.module_root / 'recovery_v4p1'

    def tearDown(self):
        self.temp.cleanup()

    def test_native_scope_is_closed_and_checker_override_forbidden(self):
        backend = native.NativeBackend325(self.frozen)
        self.assertEqual(backend.expected_table_rows, native.NATIVE_TABLE_COUNTS325)
        for config in ({'workers': 1}, {'workers': 3}, {'expected_ranks': 4},
                       {'steps': 1}, {'expected_families': 253},
                       {'arm_names': ['qual-C0', 'other']}, {'unexpected': True}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                native.NativeBackend325(self.frozen, runtime_config=config)
        with self.assertRaisesRegex(ValueError, 'native_checker_override_forbidden'):
            native.NativeBackend325(self.frozen, checker_files={'verify325.py': '/stale/copy'})

    def test_frozen_helpers_rejected_after_byte_change(self):
        package = self.root / 'frozen'
        package.mkdir()
        for name in native.FROZEN_FILES325:
            (package / name).write_bytes((self.frozen / name).read_bytes())
        (package / 'recover325.py').write_bytes(b'changed helper')
        with self.assertRaisesRegex(ValueError, 'frozen_helper_changed:recover325.py'):
            native.NativeBackend325(package)

    def test_alternate_exact_binding_copy_is_protected_input(self):
        path = self.root / 'external-binding.json'
        path.write_bytes((self.frozen / 'binding325.json').read_bytes())
        backend = native.NativeBackend325(self.frozen, binding_path=path)
        self.assertIn(path.resolve(), backend.protected_roots)

    def test_expected_provider_keys_use_both_optimized_epochs_verbatim(self):
        # Small JSON-in-gzip inputs exercise the real historical hashing reader.
        # A controlled parser supplies family dictionaries only: this is an
        # adapter choice/epoch regression, never native stream qualification.
        import gzip
        from unittest.mock import patch
        import verify325 as q

        backend = native.NativeBackend325(self.frozen)
        backend.original = self.root / 'original'
        backend.original.mkdir()
        anchors = {}
        checked = {}
        candidate = {}
        for kind, count in (('routes', 128), ('adapters', 8), ('static', 118)):
            def keys(epoch):
                return {repr(((('step', 1), ('transform.plan', epoch),
                               ('transform.scalar_adapter', epoch)), i)) +
                        ('/route/0' if kind == 'routes' else ''): {}
                        for i in range(count)}
            checked[kind] = keys(2)
            candidate[kind] = keys(1)
        for arm in native.LABELS325:
            folder = backend.original / arm
            folder.mkdir()
            path = folder / 'execution316-rank-0.bin.gz'
            path.write_bytes(gzip.compress(json.dumps(checked if arm.startswith('checked') else candidate).encode(), mtime=0))
            relative = str(path.relative_to(backend.original))
            anchors[relative] = dict(relative=relative, stored_bytes=path.stat().st_size,
                                    mtime_ns=path.stat().st_mtime_ns, sha256=native.sha256325(path),
                                    family='capture_inventory')
        calls = []
        def controlled_parser(path, rank, ranks, arm, coverage):
            calls.append((arm, rank, ranks, coverage))
            with q.v.open_stream(path) as stream:
                return json.loads(stream.read())
        original_reader = q.v.open_stream
        with patch.object(q.v, 'parse_stream', side_effect=controlled_parser):
            keys = backend._family_keys(anchors)
        self.assertIs(q.v.open_stream, original_reader)
        self.assertEqual(calls, [('qual-C0', 0, 256, 'candidate'), ('qual-L256', 0, 256, 'candidate')])
        expected = sorted((kind, key) for kind in candidate for key in candidate[kind])
        self.assertEqual(keys, expected)
        self.assertNotEqual(keys, sorted((kind, key) for kind in checked for key in checked[kind]))
        # Do not hide a genuine optimized-arm difference through normalization.
        path = backend.original / 'qual-L256/execution316-rank-0.bin.gz'
        altered = copy.deepcopy(candidate)
        old = next(iter(altered['routes']))
        altered['routes'][old + '/wrong-epoch'] = altered['routes'].pop(old)
        path.write_bytes(gzip.compress(json.dumps(altered).encode(), mtime=0))
        relative = str(path.relative_to(backend.original))
        anchors[relative].update(stored_bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                                 sha256=native.sha256325(path))
        with patch.object(q.v, 'parse_stream', side_effect=controlled_parser):
            with self.assertRaisesRegex(ValueError, 'optimized_expected_provider_key_difference'):
                backend._family_keys(anchors)
        self.assertIs(q.v.open_stream, original_reader)

    def test_every_worker_source_and_reader_attestation_is_required(self):
        backend = native.NativeBackend325(self.frozen, runtime_config={'workers': 2})
        source = {'sources': {'source': 'exact-source-hash'}}
        row = dict(passed=True, sources=source['sources'],
                   source_manifest_sha256=backend.original_binding['source_manifest_sha256'],
                   native_reader=dict(path=str(self.frozen / 'recover325.py'),
                                      sha256=native.FROZEN_FILES325['recover325.py']))
        report = {'source_bindings': {arm: [copy.deepcopy(row)] for arm in native.ARMS325}}
        self.assertEqual(backend._worker_bindings(report, source), report['source_bindings'])
        for mutation in ('missing', 'status', 'source', 'manifest', 'reader'):
            changed = copy.deepcopy(report)
            first = changed['source_bindings']['qual-C0'][0]
            if mutation == 'missing':
                changed['source_bindings']['qual-C0'] = []
            elif mutation == 'status':
                first['passed'] = False
            elif mutation == 'source':
                first['sources'] = {}
            elif mutation == 'manifest':
                first['source_manifest_sha256'] = '0' * 64
            else:
                first['native_reader']['sha256'] = '0' * 64
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                backend._worker_bindings(changed, source)

    def test_bounded_json_reserves_actual_bytes_before_open(self):
        folder = self.root / 'stage'
        folder.mkdir()
        with self.assertRaisesRegex(RuntimeError, 'phase_write_budget'):
            native.bounded_json325(folder, folder / 'proof.json', {'value': 'x' * 100}, 10)
        self.assertEqual(list(folder.iterdir()), [])
        native.bounded_json325(folder, folder / 'proof.json', {'passed': True}, 1024)
        self.assertEqual(json.loads((folder / 'proof.json').read_text()), {'passed': True})
        with self.assertRaisesRegex(ValueError, 'fresh_compact_path'):
            native.bounded_json325(folder, folder / 'proof.json', {'passed': False}, 1024)

    def test_account_writer_never_writes_outside_owned_stage_or_over_budget(self):
        folder = self.root / 'stage'
        folder.mkdir()
        module = types.SimpleNamespace(Path=Path)
        with native.bounded_account_writes325(module, folder, 10):
            module.Path(folder / 'good').write_text('small')
            with self.assertRaisesRegex(RuntimeError, 'account_write_budget'):
                module.Path(folder / 'oversize').write_text('0123456789')
            with self.assertRaisesRegex(ValueError, 'account_write_outside_owned_stage'):
                module.Path(self.root / 'protected').write_text('bad')
        self.assertIs(module.Path, Path)
        self.assertEqual((folder / 'good').read_text(), 'small')
        self.assertFalse((folder / 'oversize').exists())
        self.assertFalse((self.root / 'protected').exists())

    def test_owned_work_tree_deletion_is_narrow(self):
        folder = self.root / 'stage'
        folder.mkdir()
        other = folder / 'other'
        other.mkdir()
        (other / 'preserve').write_text('preserve')
        with self.assertRaisesRegex(ValueError, 'owned_work_tree_required'):
            native.remove_owned_tree325(folder, other)
        shadow = folder / 'work-input'
        shadow.mkdir()
        (shadow / 'rank').write_text('derived')
        native.remove_owned_tree325(folder, shadow)
        self.assertFalse(shadow.exists())
        self.assertEqual((other / 'preserve').read_text(), 'preserve')

    def test_semantic_plain_valueerror_becomes_latched_failure_class(self):
        from resume_workflow325 import SemanticFailure325, ReplayInterrupted325

        @native.semantic_phase325
        def semantic():
            raise ValueError('account325.peer_conservation')

        @native.semantic_phase325
        def interrupted():
            raise ReplayInterrupted325(15)

        @native.semantic_phase325
        def io_failure():
            raise OSError('storage unavailable')

        with self.assertRaisesRegex(SemanticFailure325, 'peer_conservation'):
            semantic()
        with self.assertRaises(ReplayInterrupted325):
            interrupted()
        with self.assertRaises(OSError):
            io_failure()

    def test_fresh_interpreter_does_not_shadow_current_modules_with_frozen_package(self):
        script = """
from pathlib import Path
import inspect,json,sys
from resume_native325 import NativeBackend325,actual_checker_paths325
assert 'parallel_replay325' not in sys.modules
root=Path.cwd()
backend=NativeBackend325(root/'recovery_v4p1')
reader=backend._activate_reader()
import parallel_replay325 as parallel
import importlib
module=importlib.import_module(reader)
assert Path(parallel.__file__).resolve()==root/'parallel_replay325.py'
assert Path(module.__file__).resolve()==root/'recovery_v4p1/recover325.py'
assert 'reader_module' in inspect.signature(parallel.cache_shards325).parameters
assert sys.path[0]!='recovery_v4p1'
print(json.dumps({'reader':reader,'paths':{k:str(v) for k,v in actual_checker_paths325().items()}}))
"""
        result = subprocess.run([sys.executable, '-B', '-c', script], cwd=self.module_root,
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['reader'], 'recovery_v4p1.recover325')

    def test_fresh_worker_hash_guard_rejects_reader_before_executing_it(self):
        package = self.root / 'suspect'
        package.mkdir()
        marker = self.root / 'executed'
        (package / 'recover325.py').write_text(f"open({str(marker)!r},'w').write('executed')\n")
        (package / 'binding325.json').write_text('{}')
        script = f"""
import sys
sys.path.append({str(self.root)!r})
from resume_native325 import native_cache_worker325
try:native_cache_worker325(('suspect.recover325',))
except ValueError as e:
 assert 'worker_reader_or_binding_bytes' in str(e),str(e)
else:raise AssertionError('unbound helper accepted')
"""
        result = subprocess.run([sys.executable, '-B', '-c', script], cwd=self.module_root,
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())

    def test_namespace_lookup_rejects_initializer_before_executing_it(self):
        package = self.root / 'suspect'
        package.mkdir()
        marker = self.root / 'initialized'
        (package / '__init__.py').write_text(f"open({str(marker)!r},'w').write('executed')\n")
        (package / 'recover325.py').write_text('unused')
        script = f"""
import sys
sys.path.append({str(self.root)!r})
from resume_native325 import frozen_reader_path325
try:frozen_reader_path325('suspect.recover325')
except ValueError as e:
 assert 'reader_namespace_required' in str(e),str(e)
else:raise AssertionError('package initializer accepted')
"""
        result = subprocess.run([sys.executable, '-B', '-c', script], cwd=self.module_root,
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())

    def test_real_pool_statuses_are_observed_and_original_factory_restored(self):
        module = types.SimpleNamespace(multiprocessing=multiprocessing)
        statuses = []
        with native.tracked_pool_context325(module, statuses):
            with module.multiprocessing.get_context('spawn').Pool(processes=2) as pool:
                self.assertEqual(pool.map(_identity325, [1, 2]), [1, 2])
        self.assertIs(module.multiprocessing, multiprocessing)
        self.assertEqual(len(statuses), 2)
        self.assertTrue(all(row['status_observed'] and row['pid'] > 0 and
                            type(row['exitcode']) is int for row in statuses))

    def test_constructor_cannot_start_any_phase_without_fresh_binding(self):
        backend = native.NativeBackend325(self.frozen)
        stage = self.root / 'stage'
        stage.mkdir()
        for method in (backend.cache, backend.work, backend.preservation):
            with self.subTest(method=method.__name__), self.assertRaisesRegex(
                    ValueError, 'fresh_binding_required'):
                method(stage, 1024**3)
        self.assertEqual(list(stage.iterdir()), [])

    def test_actual_native_work_adapter_on_source_derived_768_rank_fixture(self):
        # This exercises the real adapter and unmodified account.qualify, using
        # constructed source-derived ledgers. It is not a native qualification.
        import account325 as account
        from test_account325 import make_run
        original = self.root / 'fixture-original'
        original.mkdir()
        schema_file = self.module_root / 'work325-schema.json'
        schema = json.loads(schema_file.read_text())['records']
        (original / 'work325-schema.json').write_bytes(schema_file.read_bytes())
        maps = {arm: json.loads((self.module_root / (arm + '.json')).read_text())
                for arm in ('C0', 'L256')}
        for label in native.LABELS325:
            make_run(original / label, maps['C0' if label == 'qual-C0' else 'L256'], schema)
        initial = account.qualify(original, original)
        self.assertTrue(initial['passed'])
        ledger = {'runs': {}}
        inventory, anchors = [], {}
        for label in native.LABELS325:
            vector_path = original / label / 'work325-rank-vectors.json'
            ledger['runs'][label] = {'work_ledger': {'per_rank_artifact': {
                'bytes': vector_path.stat().st_size, 'sha256': native.sha256325(vector_path)}}}
            for rank in range(256):
                name = f'{label}/work325-rank-{rank}.tsv'
                path = original / name
                row = dict(relative=name, stored_bytes=path.stat().st_size,
                           mtime_ns=path.stat().st_mtime_ns)
                inventory.append(row)
                anchors[name] = dict(row, family='work_inventory', sha256=native.sha256325(path))
        (original / 'results.json').write_text(json.dumps(ledger))
        before = {str(path.relative_to(original)): native.sha256325(path)
                  for path in original.rglob('*') if path.is_file()}
        backend = native.NativeBackend325(self.frozen)
        # Fixture injection is confined to this test. The production workflow
        # still requires fresh binding() against the hard-pinned native records.
        backend.original = original
        backend.original_binding = dict(backend.original_binding, work_inventory=inventory,
                                        original_results=str(original))
        backend.protected_roots = (original, self.frozen)
        backend.anchors = anchors
        backend.expected_family_keys = []
        backend._last_binding = {key: '0' * 64 for key in ('source', 'checker', 'config', 'captures')}
        undersized = self.root / 'undersized-stage'
        undersized.mkdir()
        with self.assertRaisesRegex(RuntimeError, 'work_phase_preallocation'):
            backend.work(undersized, 1024**2)
        self.assertEqual(list(undersized.iterdir()), [])
        stage = self.root / 'work-stage'
        stage.mkdir()
        result = backend.work(stage, 192 * 1024**2)
        self.assertTrue(result['census']['passed'])
        self.assertEqual(result['census']['work_reads'], 768)
        self.assertFalse(result['census']['performance_qualified'])
        self.assertEqual(set(path.name for path in stage.iterdir()),
                         {'work-report.json', 'work-read-checks.json'})
        report = json.loads((stage / result['artifacts']['report']).read_text())
        receipts = json.loads((stage / result['artifacts']['actual_reads']).read_text())
        self.assertTrue(report['passed'])
        self.assertEqual(report['actual_read_hashes_matched'], 768)
        self.assertEqual(set(receipts), set(anchors))
        for name, receipt in receipts.items():
            self.assertEqual(receipt, dict(sha256=anchors[name]['sha256'],
                                           bytes=anchors[name]['stored_bytes'], complete=True))
        for arm in native.LABELS325:
            self.assertFalse(report['runs'][arm]['per_rank_artifact']['retained'])
            self.assertNotIn('path', report['runs'][arm]['per_rank_artifact'])
        self.assertNotIn(str(stage), (stage / 'work-report.json').read_text())
        after = {str(path.relative_to(original)): native.sha256325(path)
                 for path in original.rglob('*') if path.is_file()}
        self.assertEqual(before, after)


if __name__ == '__main__':
    unittest.main()
