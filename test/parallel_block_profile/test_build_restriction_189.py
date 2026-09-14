import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from build_restriction_189 import build, digest
from submit_restriction_189 import submit


class BuildRestriction189(unittest.TestCase):
    def fixture(self, root):
        kit, tree = root / 'kit', root / 'tree'
        kit.mkdir(); (tree / 'src').mkdir(parents=True)
        source = tree / 'src/example.f90'
        source.write_text('baseline\n')
        (tree / 'Makefile').write_text('build configuration\n')
        (kit / 'physics-Makefile.inc').write_text('physics include\n')
        (kit / 'stage189-candidate.patch').write_text(
            '--- a/src/example.f90\n+++ b/src/example.f90\n@@ -1 +1 @@\n-baseline\n+candidate\n')
        manifest = {'baseline_commit': 'test-base', 'baseline_sources': {'src/example.f90': digest(source)},
                    'candidate_sources': {'src/example.f90': hashlib.sha256(b'candidate\n').hexdigest()},
                    'makefile_sha256': digest(tree / 'Makefile'),
                    'candidate_patch_sha256': digest(kit / 'stage189-candidate.patch')}
        (kit / 'stage189-sources.json').write_text(json.dumps(manifest))
        return kit, tree

    def fake_run(self, tree, mutate_physics=False):
        real_run = subprocess.run
        def run(command, **kwargs):
            if command[0] == 'make':
                self.assertIn('FP_CONTRACT=off', command)
                self.assertIn('PARAM=param_J5', command)
                self.assertIn('-j1', command)
                name = next(c.split('=', 1)[1] for c in command if c.startswith('BIN_DIR='))
                expected = 'baseline\n' if name.endswith('baseline') else 'candidate\n'
                self.assertEqual((tree / 'src/example.f90').read_text(), expected)
                (tree / name).mkdir()
                (tree / name / 'climate').write_text(name)
                objects = tree / 'src/physics/simple_physics/phyparam/obj'
                objects.mkdir(parents=True, exist_ok=True)
                for i in range(23):
                    (objects / f'{i}.o').write_text(str(i))
                if mutate_physics and name.endswith('candidate'):
                    (objects / '0.o').write_text('changed')
                return subprocess.CompletedProcess(command, 0)
            if command[0] in ('bash', 'mpirun'):
                return subprocess.CompletedProcess(command, 0, stdout='MPI version', stderr='')
            return real_run(command, **kwargs)
        return run

    def test_build_order_flags_and_physics_guard(self):
        for changed in (False, True):
            with self.subTest(changed=changed), tempfile.TemporaryDirectory() as tmp:
                kit, tree = self.fixture(Path(tmp))
                with patch('build_restriction_189.subprocess.run', side_effect=self.fake_run(tree, changed)), \
                        patch('build_restriction_189.subprocess.check_output', return_value='test compiler'):
                    if changed:
                        with self.assertRaisesRegex(ValueError, 'identical external physics'):
                            build(kit, tree)
                    else:
                        build(kit, tree)
                record = json.loads((tree / 'restriction-189-builds.json').read_text())
                self.assertEqual(record['status'], 'failed' if changed else 'passed')
                if not changed:
                    self.assertEqual(list(record['builds']), ['baseline', 'candidate', 'checked'])
                    self.assertEqual(len(record['physics_objects']), 23)

    def test_stale_build_is_rejected_before_any_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit, tree = self.fixture(Path(tmp))
            (tree / 'build-189-baseline').mkdir()
            with patch('build_restriction_189.subprocess.run') as run:
                with self.assertRaisesRegex(ValueError, 'Fresh build'):
                    build(kit, tree)
                run.assert_not_called()

    def test_submission_checks_binary_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); kit, tree = self.fixture(root)
            with patch('build_restriction_189.subprocess.run', side_effect=self.fake_run(tree)), \
                    patch('build_restriction_189.subprocess.check_output', return_value='test compiler'):
                build(kit, tree)
            fixture = root / 'fixture'; fixture.mkdir()
            with patch('submit_restriction_189.subprocess.run') as run:
                submit(kit, tree, fixture, root / 'out')
                self.assertIn('--profiles', run.call_args.args[0])
            (tree / 'bin-189-candidate/climate').write_text('stale replacement')
            with patch('submit_restriction_189.subprocess.run') as run:
                with self.assertRaisesRegex(ValueError, 'Executable changed'):
                    submit(kit, tree, fixture, root / 'out')
                run.assert_not_called()
