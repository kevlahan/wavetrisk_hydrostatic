"""Offline single-submission, immutable-package and bounded-collection tests."""
import hashlib
import io
import json
import os
import signal
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time
import types
import unittest
from unittest.mock import patch

import control_native325 as control
import build_native325 as builder


class NativeCampaign325Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()

    def tearDown(self):
        self.temp.cleanup()

    def kit(self, path=None):
        path = self.root / 'kit' if path is None else path
        path.mkdir()
        (path / 'helper.py').write_text('verified source\n')
        manifest = {'helper.py': control.sha(path / 'helper.py')}
        (path / 'native-files-sha256.json').write_text(json.dumps(manifest))
        return path

    def state(self):
        kit = self.kit()
        work = self.root / 'work'
        work.mkdir()
        return dict(status='submitted', job=999, kit=str(kit), work=str(work), out=str(work / 'out'),
                    package_manifest_sha256=control.package_check(kit))

    def test_package_inventory_and_same_size_mutation(self):
        kit = self.kit()
        before = control.package_check(kit)
        self.assertEqual(len(before), 64)
        p = kit / 'helper.py'
        p.write_bytes(b'x' * p.stat().st_size)
        with self.assertRaisesRegex(ValueError, 'package_bytes'):
            control.package_check(kit)

    def test_unknown_file_and_symlink_fail_package_gate(self):
        kit = self.kit()
        (kit / 'unexpected').write_text('extra')
        with self.assertRaisesRegex(ValueError, 'package_inventory'):
            control.package_check(kit)
        (kit / 'unexpected').unlink()
        (kit / 'link').symlink_to(kit / 'helper.py')
        with self.assertRaisesRegex(ValueError, 'package_entry_type'):
            control.package_check(kit)

    def test_unlisted_bytecode_is_not_ignored(self):
        kit = self.kit()
        bytecode = kit / '__pycache__'
        bytecode.mkdir()
        (bytecode / 'injected.pyc').write_bytes(b'not sealed')
        with self.assertRaisesRegex(ValueError, 'package_inventory'):
            control.package_check(kit)

    def test_expired_squeue_job_falls_back_to_terminal_accounting(self):
        answers = [types.SimpleNamespace(returncode=1, stdout='', stderr='Invalid job id'),
                   types.SimpleNamespace(returncode=0, stdout='999|COMPLETED|\n', stderr='')]
        with patch.object(control.subprocess, 'run', side_effect=answers):
            self.assertFalse(control.queued(999))

    def test_scheduler_errors_and_unknown_accounting_fail_closed(self):
        for answer in (types.SimpleNamespace(returncode=1, stdout='', stderr='offline'),
                       types.SimpleNamespace(returncode=0, stdout='', stderr='')):
            with self.subTest(answer=answer), patch.object(control.subprocess, 'run', return_value=answer):
                with self.assertRaisesRegex(RuntimeError, 'Scheduler state unknown'):
                    control.queued(999)
        with self.assertRaises(ValueError):
            control.queued(True)

    def test_active_accounting_is_not_collectable(self):
        answers = [types.SimpleNamespace(returncode=0, stdout='', stderr=''),
                   types.SimpleNamespace(returncode=0, stdout='999|RUNNING|\n', stderr='')]
        with patch.object(control.subprocess, 'run', side_effect=answers):
            self.assertTrue(control.queued(999))
        state = self.state()
        with patch.object(control, 'queued', return_value=True), self.assertRaisesRegex(
                ValueError, 'terminal_job_required'):
            control.collect(state)

    def test_collection_omits_large_native_data_without_reading_it(self):
        state = self.state()
        out = Path(state['out'])
        (out / 'state/phases/cache').mkdir(parents=True)
        (out / 'state/partial/work/work-input/qual-C0').mkdir(parents=True)
        protected = [out / 'state/phases/cache/left.sqlite',
                     out / 'state/partial/work/work-input/qual-C0/work325-rank-0.tsv',
                     out / 'state/partial/work/work-input/qual-C0/work325-rank-vectors.json']
        for path in protected:
            with path.open('wb') as stream:
                stream.truncate(4 * 1024**3)
        (out / 'supervisor-report.json').write_text('{"status":"failed"}')
        original = Path.read_bytes
        def guarded(path):
            if path in protected:
                raise AssertionError('Excluded native content was read')
            return original(path)
        with patch.object(Path, 'read_bytes', guarded):
            files, excluded, count = control.selected_files(state)
        self.assertIn('work/out/supervisor-report.json', [name for _, name, _, _ in files])
        self.assertEqual(len(excluded), 3)
        self.assertLess(count, 1024)

    def test_collection_roundtrip_and_identity_reuse(self):
        state = self.state()
        out = Path(state['out'])
        (out / 'state/phases/preservation').mkdir(parents=True)
        (out / 'state/phases/preservation/preservation.json').write_text('{"passed":true}')
        accounting = types.SimpleNamespace(returncode=0, stdout='999|COMPLETED|0:0|start|end|00:01:00|60|bb13|\n',
                                          stderr='', args=['sacct'])
        with patch.object(control.Path, 'home', return_value=self.root), patch.object(
                control, 'queued', return_value=False), patch.object(control.subprocess, 'run', return_value=accounting):
            archive = control.collect(state)
            self.assertEqual(control.collect(state), archive)
            self.assertGreater(control.verify_archive(archive), 5)
            changed = dict(state, package_manifest_sha256='0' * 64)
            with self.assertRaisesRegex(ValueError, 'existing_review_identity'):
                control.collect(changed)
        with tarfile.open(archive) as contents:
            self.assertTrue(all(member.isfile() for member in contents.getmembers()))
            self.assertIn('review/work/out/state/phases/preservation/preservation.json', contents.getnames())

    def test_no_duplicate_submission_when_state_exists(self):
        kit = self.kit()
        latest = self.root / (control.PREFIX + '-latest.json')
        latest.write_text('{}')
        with patch.object(control.Path, 'home', return_value=self.root), patch.object(
                control, 'queued', return_value=False), patch.object(control.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'state_exists_no_automatic_retry'):
                control.submit(kit)
        run.assert_not_called()

    def test_failed_v5_must_be_confirmed_terminal_before_new_submission(self):
        kit = self.kit()
        with patch.object(control.Path, 'home', return_value=self.root), patch.object(
                control, 'queued', side_effect=lambda job: job == 55575) as queued, patch.object(
                control.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'predecessor_still_active:55575'):
                control.submit(kit)
        self.assertEqual([call.args[0] for call in queued.call_args_list],
                         [55569, 55571, 55572, 55573, 55574, 55575])
        run.assert_not_called()
        self.assertFalse((self.root / (control.PREFIX + '-latest.json')).exists())

    def test_preflight_failure_blocks_sbatch_and_retains_state(self):
        kit = self.kit()
        with patch.object(control.Path, 'home', return_value=self.root), patch.object(
                control, 'queued', return_value=False), patch.object(control, 'run_preflight',
                side_effect=ValueError('unit tests failed')), patch.object(control.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'unit tests failed'):
                control.submit(kit)
        run.assert_not_called()
        state = json.loads((self.root / (control.PREFIX + '-latest.json')).read_text())
        self.assertEqual(state['status'], 'preflight_failed')
        self.assertIsNone(state['job'])

    def test_single_verified_submission_has_private_work_and_no_resume_flags(self):
        kit = self.kit()
        response = types.SimpleNamespace(returncode=0, stdout='999;cluster\n', stderr='')
        with patch.object(control.Path, 'home', return_value=self.root), patch.object(
                control, 'queued', return_value=False), patch.object(control, 'run_preflight',
                return_value={'passed': True}), patch.object(control.subprocess, 'run', return_value=response) as run:
            state = control.submit(kit)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.args[0][0], 'sbatch')
        self.assertEqual(state['job'], 999)
        self.assertEqual(state['predecessor_job'], 55575)
        self.assertEqual(state['predecessor_elapsed_seconds'], 0)
        self.assertEqual(state['new_budget_remaining_seconds'], 1800)
        self.assertEqual(Path(state['work']).stat().st_mode & 0o777, 0o700)
        self.assertFalse(Path(state['out']).exists())
        self.assertIn('--output=' + str(Path(state['work']) / 'slurm.log'), run.call_args.args[0])

    def test_archive_path_traversal_rejected(self):
        path = self.root / 'unsafe.tar.gz'
        with tarfile.open(path, 'w:gz') as archive:
            info = tarfile.TarInfo('../escape')
            info.size = 1
            archive.addfile(info, io.BytesIO(b'x'))
        with self.assertRaisesRegex(ValueError, 'review_member_safety'):
            control.verify_archive(path)

    def test_frozen_payload_requires_exact_historical_manifest(self):
        kit = self.kit()
        (kit / 'native-files-sha256.json').rename(kit / 'files-sha256.json')
        with self.assertRaisesRegex(ValueError, 'historical_manifest_identity'):
            builder.frozen_payload(self.root, kit)
        manifest_sha = control.sha(kit / 'files-sha256.json')
        with patch.object(builder, 'FROZEN_MANIFEST_SHA', manifest_sha):
            payload = builder.frozen_payload(self.root, kit)
        self.assertEqual(set(payload), {'helper.py', 'files-sha256.json'})

    def test_submission_script_is_one_bound_python_verifier_step(self):
        # Source tree and sealed package place the script beside this test.
        script = (Path(__file__).parent / 'verification.sbatch').read_text()
        self.assertIn('#SBATCH --time=00:30:00', script)
        self.assertIn('#SBATCH --exclusive', script)
        self.assertIn('--cpu-bind=cores python3 -B', script)
        self.assertIn('native_runner325.py" supervise', script)
        self.assertNotIn('#SBATCH --signal', script)
        self.assertNotIn('mpirun', script)

    def test_preflight_bounds_log_before_write(self):
        kit = self.root / 'preflight-kit'
        kit.mkdir()
        (kit / 'test_small.py').write_text('import unittest\nclass Small(unittest.TestCase):\n def test_pass(self):pass\n')
        log = self.root / 'test.log'
        result = control.run_preflight(kit, log)
        self.assertTrue(result['passed'])
        self.assertLess(log.stat().st_size, control.MAX_TEST_LOG)
        (kit / 'test_small.py').write_text('import unittest\nclass Small(unittest.TestCase):\n def test_pass(self):print("x"*1024)\n')
        with patch.object(control, 'MAX_TEST_LOG', 64), self.assertRaisesRegex(ValueError, 'preflight_log_cap'):
            control.run_preflight(kit, self.root / 'limited.log')
        self.assertLessEqual((self.root / 'limited.log').stat().st_size, 64)

    def test_preflight_kills_ignoring_grandchild_after_parent_exits(self):
        kit = self.root / 'orphan-kit'
        kit.mkdir()
        pid_file = self.root / 'orphan.pid'
        heartbeat = self.root / 'orphan-heartbeat'
        child = ("import os,signal,time\nsignal.signal(signal.SIGTERM,signal.SIG_IGN)\n"
                 f"open({str(pid_file)!r},'w').write(str(os.getpid()))\n"
                 "end=time.monotonic()+60\nwhile time.monotonic()<end:\n"
                 f" open({str(heartbeat)!r},'w').write(str(time.time_ns()))\n time.sleep(.01)\n")
        source = ('import subprocess,sys,time,unittest\nfrom pathlib import Path\n'
                  'class Orphan(unittest.TestCase):\n def test_parent_exits(self):\n'
                  f'  subprocess.Popen([sys.executable,"-c",{child!r}])\n'
                  f'  path=Path({str(pid_file)!r})\n'
                  '  deadline=time.monotonic()+2\n'
                  '  while not path.exists() and time.monotonic()<deadline:time.sleep(.01)\n'
                  '  self.assertTrue(path.exists())\n')
        (kit / 'test_orphan.py').write_text(source)
        pid = None
        try:
            with self.assertRaisesRegex(ValueError, 'preflight_timeout'):
                control.run_preflight(kit, self.root / 'orphan.log', timeout=.5)
            pid = int(pid_file.read_text())
            # Check actual execution rather than signal0/ps enumeration, which
            # managed local sandboxes can forbid even for an owned descendant.
            before = heartbeat.read_bytes()
            time.sleep(.1)
            self.assertEqual(heartbeat.read_bytes(), before)
        finally:
            if pid is None and pid_file.exists():
                pid = int(pid_file.read_text())
            if pid is not None:
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass


if __name__ == '__main__':
    unittest.main()
