"""Adversarial byte/certificate and real process-interruption store tests."""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest

from resume_state325 import PhaseStore325, ResumeStateError325, PHASES325


BINDING = {k: hashlib.sha256(k.encode()).hexdigest()
           for k in ('source', 'checker', 'config', 'captures')}


class ResumeState325Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name).resolve()
        self.root = self.base / 'derived'
        self.inputs = self.base / 'inputs'
        self.inputs.mkdir()
        (self.inputs / 'capture.gz').write_bytes(b'protected-original')
        self.options = dict(protected_roots=[self.inputs], metadata_limit=4096,
                            cap_bytes=1024 * 1024)

    def tearDown(self):
        self.temp.cleanup()

    def store(self, *, resume=False, binding=None):
        return PhaseStore325(self.root, BINDING if binding is None else binding,
                             explicit_resume=resume, **self.options)

    def complete(self, store, phase, value=b'complete-artifact'):
        stage = store.start(phase)
        store.reserve(phase, len(value) + 100)
        (stage / 'proof.bin').write_bytes(value)
        return store.commit(phase, {'proof': 'proof.bin'}, census={'items': 1},
                            notes=['closed complete artifact'])

    def child(self, source, *, wait=True):
        header = ("from pathlib import Path\nimport os, signal, time\n"
                  "from resume_state325 import PhaseStore325\n"
                  f"binding={BINDING!r}\nroot=Path({str(self.root)!r})\n"
                  f"options={dict(self.options, protected_roots=[str(self.inputs)])!r}\n")
        process = subprocess.Popen([sys.executable, '-c', header + source],
                                   cwd=Path(__file__).resolve().parent,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if wait:
            output, error = process.communicate(timeout=10)
            return process.returncode, output, error
        return process

    def test_full_dag_resume_rehash_and_root_budget(self):
        with self.store() as store:
            initial = store.available_bytes()
            for phase in PHASES325:
                result = self.complete(store, phase)
                self.assertEqual(result['certificate']['phase'], phase)
                self.assertEqual(store.artifact(phase, 'proof').read_bytes(), b'complete-artifact')
            self.assertEqual(initial - store.available_bytes(), 6 * len(b'complete-artifact'))
        with self.store(resume=True) as resumed:
            self.assertEqual(set(resumed.validate_all()), set(PHASES325))

    def test_explicit_resume_required(self):
        with self.store():
            pass
        with self.assertRaisesRegex(ResumeStateError325, 'explicit_resume_required'):
            with self.store():
                pass

    def test_every_binding_digest_invalidates_state(self):
        with self.store() as store:
            self.complete(store, 'cache')
        for key in BINDING:
            changed = dict(BINDING, **{key: '0' * 64})
            with self.subTest(key=key), self.assertRaisesRegex(
                    ResumeStateError325, 'identity_or_binding_changed'):
                with self.store(resume=True, binding=changed):
                    pass

    def test_same_size_and_mtime_artifact_mutation_rejected(self):
        with self.store() as store:
            self.complete(store, 'cache')
            path = store.artifact('cache', 'proof')
        before = path.stat()
        path.write_bytes(b'X' * before.st_size)
        os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
        with self.assertRaisesRegex(ResumeStateError325, 'artifact_bytes_changed'):
            with self.store(resume=True):
                pass

    def test_certificate_exact_fields_and_dependency_digest(self):
        with self.store() as store:
            self.complete(store, 'cache')
            self.complete(store, 'providers-C0')
        cert = self.root / 'phases' / 'cache' / 'certificate.json'
        original = cert.read_bytes()
        obj = json.loads(original)
        obj['ignore_this'] = True
        cert.write_text(json.dumps(obj))
        with self.assertRaisesRegex(ResumeStateError325, 'certificate_fields'):
            with self.store(resume=True):
                pass
        cert.write_bytes(original)
        obj = json.loads(original)
        obj['notes'].append('changed dependency certificate')
        cert.write_text(json.dumps(obj))
        with self.assertRaisesRegex(ResumeStateError325, 'dependency_certificate_changed'):
            with self.store(resume=True):
                pass

    def test_failure_latch_blocks_current_and_future_session(self):
        with self.store() as store:
            self.complete(store, 'cache')
            store.mark_semantic_failed('provider mismatch')
            with self.assertRaisesRegex(ResumeStateError325, 'semantic_failed_latch'):
                store.start('providers-C0')
            with self.assertRaisesRegex(ResumeStateError325, 'semantic_failed_latch'):
                store.mark_interrupted(signal.SIGTERM)
        with self.assertRaisesRegex(ResumeStateError325, 'semantic_failed_latch'):
            with self.store(resume=True):
                pass

    def test_partial_never_counts_as_complete_and_explicit_discard(self):
        with self.store() as store:
            stage = store.start('cache')
            store.reserve('cache', 100)
            (stage / 'partial.bin').write_bytes(b'partial')
            store.mark_interrupted(signal.SIGTERM)
        with self.store(resume=True) as store:
            self.assertFalse(store.completed('cache'))
            with self.assertRaisesRegex(ResumeStateError325, 'explicit_partial_discard_required'):
                store.start('cache')
            stage = store.start('cache', discard_partial=True)
            self.assertEqual(list(stage.iterdir()), [])
            store.reserve('cache', 100)
            (stage / 'proof.bin').write_bytes(b'new')
            store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
        self.assertEqual((self.inputs / 'capture.gz').read_bytes(), b'protected-original')

    def test_unknown_phase_and_missing_dependency_fail(self):
        with self.store() as store:
            with self.assertRaisesRegex(ResumeStateError325, 'unknown_phase'):
                store.start('fake-phase')
            with self.assertRaisesRegex(ResumeStateError325, 'missing_dependency'):
                store.start('pair')

    def test_progress_record_and_partial_certificate_never_prove_completion(self):
        with self.store():
            pass
        path = self.root / 'session.json'
        obj = json.loads(path.read_text())
        obj.update(status='complete', phase='preservation')
        path.write_text(json.dumps(obj))
        with self.store(resume=True) as store:
            self.assertFalse(store.completed('cache'))
            self.assertFalse(store.completed('preservation'))
            with self.assertRaisesRegex(ResumeStateError325, 'missing_dependency'):
                store.start('pair')

    def test_database_corruption_rehashed_even_when_census_and_stat_unchanged(self):
        import sqlite3
        with self.store() as store:
            stage = store.start('cache')
            store.reserve('cache', 65536)
            database = stage / 'cache.sqlite'
            connection = sqlite3.connect(database)
            connection.execute('create table proof(value text)')
            connection.execute("insert into proof values ('retained')")
            connection.commit()
            connection.close()
            store.commit('cache', {'database': 'cache.sqlite'}, census={'proof_rows': 1})
            database = store.artifact('cache', 'database')
        before = database.stat()
        data = database.read_bytes()
        self.assertIn(b'retained', data)
        database.write_bytes(data.replace(b'retained', b'corrupt!', 1))
        os.utime(database, ns=(before.st_atime_ns, before.st_mtime_ns))
        self.assertEqual(database.stat().st_size, before.st_size)
        with self.assertRaisesRegex(ResumeStateError325, 'artifact_bytes_changed'):
            with self.store(resume=True):
                pass

    def test_extra_file_and_unsafe_artifact_paths_fail(self):
        with self.store() as store:
            stage = store.start('cache')
            store.reserve('cache', 100)
            (stage / 'proof.bin').write_bytes(b'proof')
            (stage / 'database-journal').write_bytes(b'leftover')
            with self.assertRaisesRegex(ResumeStateError325, 'artifact_census'):
                store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
            for value in ('../inputs/capture.gz', '/absolute', 'nested/../proof.bin',
                          'certificate.json', './proof.bin'):
                with self.subTest(path=value), self.assertRaisesRegex(
                        ResumeStateError325, '(relative_path|reserved_artifact_path)'):
                    store.commit('cache', {'proof': value}, census={'items': 1})

    def test_leaf_parent_symlink_and_hardlink_refused(self):
        alias = self.base / 'alias'
        alias.symlink_to(self.inputs, target_is_directory=True)
        with self.assertRaisesRegex(ResumeStateError325, 'symlink_path'):
            PhaseStore325(alias / 'derived', BINDING, **self.options)
        with self.store() as store:
            stage = store.start('cache')
            store.reserve('cache', 100)
            bad = stage / 'proof.bin'
            bad.symlink_to(self.inputs / 'capture.gz')
            with self.assertRaisesRegex(ResumeStateError325, 'symlink_in_store'):
                store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
            bad.unlink()
            os.link(self.inputs / 'capture.gz', bad)
            with self.assertRaisesRegex(ResumeStateError325, 'hardlink_in_store'):
                store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})

    def test_protected_overlap_and_path_traversal(self):
        for root in (self.inputs / 'new', self.base):
            with self.subTest(root=root), self.assertRaisesRegex(
                    ResumeStateError325, 'protected_root_overlap'):
                PhaseStore325(root, BINDING, **self.options)
        with self.assertRaisesRegex(ResumeStateError325, 'absolute_clean_path'):
            PhaseStore325(self.base / 'x' / '..' / 'derived', BINDING, **self.options)

    def test_reserve_before_write_and_all_retained_bytes_accounted(self):
        with self.store() as store:
            stage = store.start('cache')
            with self.assertRaisesRegex(ResumeStateError325, 'reservation_cap'):
                store.reserve('cache', store.available_bytes() + 1)
            store.reserve('cache', 200)
            (stage / 'proof.bin').write_bytes(b'x' * 201)
            with self.assertRaisesRegex(ResumeStateError325, 'unreserved_artifact_bytes'):
                store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
            (stage / 'proof.bin').write_bytes(b'x' * 200)
            store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
            self.assertEqual(store.available_bytes(), store.payload_cap - 200)
            stage = store.start('providers-C0')
            with self.assertRaisesRegex(ResumeStateError325, 'reservation_cap'):
                store.reserve('providers-C0', store.payload_cap - 199)
            self.assertFalse((stage / 'should-not-exist').exists())

    def test_unreserved_partial_bytes_rejected_on_resume(self):
        with self.store() as store:
            stage = store.start('cache')
            (stage / 'unreserved.bin').write_bytes(b'x')
        with self.assertRaisesRegex(ResumeStateError325, 'partial_allocation_exceeded'):
            with self.store(resume=True):
                pass

    def test_metadata_schema_duplicate_unknown_path_and_bound(self):
        with self.store():
            pass
        session = self.root / 'session.json'
        original = session.read_bytes()
        session.write_text('{"schema":1,"schema":1}')
        with self.assertRaisesRegex(ResumeStateError325, 'duplicate_json_key'):
            with self.store(resume=True):
                pass
        session.write_bytes(original)
        (self.root / 'phases' / 'unexpected').mkdir()
        with self.assertRaisesRegex(ResumeStateError325, 'unknown_phase_path'):
            with self.store(resume=True):
                pass
        (self.root / 'phases' / 'unexpected').rmdir()
        session.write_bytes(b' ' * 4097)
        with self.assertRaisesRegex(ResumeStateError325, 'metadata_bound'):
            with self.store(resume=True):
                pass

    def test_real_sigkill_lock_released_and_last_certificate_reusable(self):
        process = self.child("""
with PhaseStore325(root,binding,**options) as store:
    stage=store.start('cache'); store.reserve('cache',100)
    (stage/'proof.bin').write_bytes(b'completed')
    store.commit('cache',{'proof':'proof.bin'},census={'items':1})
    print('locked',flush=True)
    time.sleep(60)
""", wait=False)
        try:
            self.assertEqual(process.stdout.readline().strip(), 'locked')
            with self.assertRaisesRegex(ResumeStateError325, 'concurrent_writer'):
                with self.store(resume=True):
                    pass
            process.kill()
            process.communicate(timeout=10)
            self.assertEqual(process.returncode, -signal.SIGKILL)
            with self.store(resume=True) as store:
                self.assertTrue(store.completed('cache'))
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()

    def test_real_sigkill_at_certificate_and_directory_commit_boundaries(self):
        # A fsynced temporary certificate is never authority. Neither is a
        # complete certificate still inside the uncommitted staging directory.
        # The directory rename atomically publishes both artifact and certificate.
        for boundary, complete in (('certificate-before-replace', False),
                                   ('directory-before-rename', False),
                                   ('directory-after-rename', True)):
            with self.subTest(boundary=boundary):
                self.root = self.base / boundary
                source = f"""
import resume_state325 as module
boundary={boundary!r}
original_replace=module.os.replace
original_rename=module.os.rename
def replace(source,destination):
    if boundary=='certificate-before-replace' and Path(destination).name=='certificate.json':
        os.kill(os.getpid(),signal.SIGKILL)
    return original_replace(source,destination)
def rename(source,destination):
    if boundary=='directory-before-rename': os.kill(os.getpid(),signal.SIGKILL)
    value=original_rename(source,destination)
    if boundary=='directory-after-rename': os.kill(os.getpid(),signal.SIGKILL)
    return value
module.os.replace=replace
module.os.rename=rename
with PhaseStore325(root,binding,**options) as store:
    stage=store.start('cache'); store.reserve('cache',100)
    (stage/'proof.bin').write_bytes(b'completed')
    store.commit('cache',{{'proof':'proof.bin'}},census={{'items':1}})
"""
                code, output, error = self.child(source)
                self.assertEqual(code, -signal.SIGKILL, (output, error))
                with self.store(resume=True) as store:
                    self.assertEqual(store.completed('cache'), complete)
                    if not complete:
                        with self.assertRaisesRegex(ResumeStateError325,
                                                    'explicit_partial_discard_required'):
                            store.start('cache')
                        stage = store.start('cache', discard_partial=True)
                        self.assertEqual(list(stage.iterdir()), [])
                        store.reserve('cache', 100)
                        (stage / 'proof.bin').write_bytes(b'restarted')
                        store.commit('cache', {'proof': 'proof.bin'}, census={'items': 1})
                    self.assertEqual(set(store.validate_all()), {'cache'})

    def test_real_sigkill_during_semantic_failure_recording_still_blocks_resume(self):
        source = """
import resume_state325 as module
original_fsync_dir=module._fsync_dir
def fsync_dir(path):
    value=original_fsync_dir(path)
    if (root/'semantic-failure.json').exists(): os.kill(os.getpid(),signal.SIGKILL)
    return value
with PhaseStore325(root,binding,**options) as store:
    stage=store.start('cache'); store.reserve('cache',100)
    (stage/'proof.bin').write_bytes(b'complete')
    store.commit('cache',{'proof':'proof.bin'},census={'items':1})
    module._fsync_dir=fsync_dir
    store.mark_semantic_failed('bad provider')
"""
        code, output, error = self.child(source)
        self.assertEqual(code, -signal.SIGKILL, (output, error))
        self.assertTrue((self.root / 'semantic-failure.json').exists())
        with self.assertRaisesRegex(ResumeStateError325, 'metadata_json'):
            with self.store(resume=True):
                pass


if __name__ == '__main__':
    unittest.main()
