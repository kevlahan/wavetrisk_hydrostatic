"""Exercise cluster orchestration locally with a fake launcher, never MPI."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import cluster_187 as campaign

LOG = '''min_level = 5 DOMAIN_LEVEL = 2
Restarting from checkpoint 3
00000322 0.3100 d dt = 207.1 s Jmax = 7 dof = 999 cpu = 2.0E+00
Saving checkpoint 4
Restarting from checkpoint 4
00000322 0.3210 d dt = 207.1 s Jmax = 7 dof = 999 cpu = 3.0E+00
Total cpu time = 5.0
'''


def fixture(root):
    root.mkdir()
    (root / 'grids').mkdir()
    (root / 'grids' / 'HR95JT_WT_004').write_text('immutable grid')
    (root / 'simple.in').write_text('''run_id simple
physics_type Simple
max_level 7
zlevels 30
Nsoil 10
NCAR_topo F
sso F
resume 3
CP_EVERY 1
time_end 0.3225
dt_write 0.01
''')
    n = 160
    payload = struct.pack('<idqi', 322, 0.3, 10000, 3)
    payload += np.zeros(3 * 41 + 41 * 80, dtype='<f8').tobytes()
    sizes = np.full(n, len(payload), dtype='<i8')
    offsets = np.arange(n, dtype='<i8') * len(payload) + 24 + 20 * n
    raw = struct.pack('<3q', 0x5741564554524953, 1, n)
    raw += np.zeros(n, dtype='<i4').tobytes() + offsets.tobytes() + sizes.tobytes() + payload * n
    subprocess.run(['zstd', '-q', '-o', str(root / 'simple_checkpoint_0003.bin.zst')], input=raw, check=True)


class ClusterProtocolTests(unittest.TestCase):
    def test_schedule_excludes_warmups_and_profiles(self):
        jobs = campaign.schedule(3, True)
        self.assertEqual([x[1] for x in jobs if x[4] == 'timing'],
                         ['baseline', 'candidate', 'candidate', 'baseline', 'baseline', 'candidate'])
        self.assertEqual(sum(x[4] == 'warmup' for x in jobs), 2)
        self.assertEqual(sum(x[4] == 'profile' for x in jobs), 2)
        self.assertEqual([x[2] for x in jobs[:2]], [False, True])

    def test_coverage_rejects_missing_restart_step_extra_checkpoint_and_wrong_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for index in (3, 4):
                (root / f'simple_checkpoint_{index:04d}.bin.zst').touch()
            log = root / 'run.log'
            log.write_text(LOG)
            _, evidence = campaign.validate_log(root, 'simple', 3, [4])
            self.assertFalse(evidence['post_restart_remaps']['4'])
            with self.assertRaisesRegex(ValueError, 'No remap'):
                campaign.validate_log(root, 'simple', 3, [4], True)
            log.write_text(LOG.replace('00000322 0.3210 d dt = 207.1 s Jmax = 7 dof = 999 cpu = 3.0E+00\n', ''))
            with self.assertRaisesRegex(ValueError, 'No timestep'):
                campaign.validate_log(root, 'simple', 3, [4])
            log.write_text(LOG.replace('min_level = 5', 'min_level = 4'))
            with self.assertRaisesRegex(ValueError, 'param_J5'):
                campaign.validate_log(root, 'simple', 3, [4])
            log.write_text(LOG)
            (root / 'simple_checkpoint_0005.bin.zst').touch()
            with self.assertRaisesRegex(ValueError, 'checkpoint set'):
                campaign.validate_log(root, 'simple', 3, [4])

    @unittest.skipUnless(shutil.which('zstd'), 'zstd required')
    def test_complete_campaign_and_fail_closed_runtime(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root / 'fixture'
            fixture(seed)
            original = campaign.digest(seed / 'simple_checkpoint_0003.bin.zst')
            launch = root / 'srun'
            launch.write_text('#!/usr/bin/env python3\nimport os,sys\ni=sys.argv.index("./climateJ5")\nos.execv(sys.argv[i],sys.argv[i:])\n')
            launch.chmod(0o755)
            binaries = {}
            for kind in ('baseline', 'candidate', 'checked'):
                binary = root / kind
                binary.write_text('#!/usr/bin/env python3\nfrom pathlib import Path\nimport shutil,os\n'
                                  'shutil.copy2("simple_checkpoint_0003.bin.zst","simple_checkpoint_0004.bin.zst")\n'
                                  'shutil.copy2("simple_checkpoint_0003.bin.zst","simple_checkpoint_0005.bin.zst")\n'
                                  'Path("env.json").write_text(__import__("json").dumps(dict(os.environ)))\n'
                                  'assert not Path("grids").is_symlink()\n'
                                  'assert Path("grids/HR95JT_WT_004").read_text()=="immutable grid"\n'
                                  'if Path.cwd().name == "checked-off":\n'
                                  ' shutil.move(' + repr(str(seed / 'grids')) + ',' + repr(str(root / 'original-grids-moved')) + ')\n'
                                  'print(' + repr(LOG.replace('Total cpu time = 5.0', 'Saving checkpoint 5\nRestarting from checkpoint 5\n00000322 0.3212 d dt = 103.6 s Jmax = 7 dof = 999 cpu = 4.0E+00\nTotal cpu time = 9.0')) + ')\n')
                binary.chmod(0o755)
                binaries[kind] = binary
            args = argparse.Namespace(**binaries, fixture=seed, out=root / 'results', ranks=83,
                                      repeats=1, profiles=True, expected_checkpoints=[4, 5],
                                      require_post_restart_remap=False, build_notes='Synthetic harness test only')
            env = {'PATH': str(root) + os.pathsep + os.environ['PATH'], 'SLURM_JOB_ID': 'test',
                   'SLURM_NTASKS': '83', 'WAVETRISK_UNKNOWN_PROBE': '1', 'OMP_NUM_THREADS': '8'}
            with patch.dict(os.environ, env), contextlib.redirect_stdout(io.StringIO()):
                campaign.run(args)
            result = json.loads((args.out / 'results.json').read_text())
            self.assertEqual(result['status'], 'passed')
            self.assertEqual(len(result['runs']), 8)
            self.assertEqual(len(result['summary']['candidate']['elapsed_seconds']), 1)
            self.assertTrue(all(g['passed'] for g in result['gates'].values()))
            actual_env = json.loads((args.out / 'checked-on' / 'env.json').read_text())
            self.assertNotIn('WAVETRISK_UNKNOWN_PROBE', actual_env)
            self.assertEqual(actual_env['OMP_NUM_THREADS'], '1')
            self.assertEqual(actual_env['WAVETRISK_VALIDATE_BLOCK_DYNAMICS'], '1')
            self.assertEqual(campaign.digest(seed / 'simple_checkpoint_0003.bin.zst'), original)
            # The later runs used copies even after the original grids moved.
            self.assertFalse((seed / 'grids').exists())
            shutil.move(root / 'original-grids-moved', seed / 'grids')
            # A zero-exit runtime error must stop before timing starts.
            binaries['checked'].write_text('#!/usr/bin/env python3\nprint("ERROR STOP injected")\n')
            args.out = root / 'failed'
            with patch.dict(os.environ, env), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'Incomplete'):
                    campaign.run(args)
            failed = json.loads((args.out / 'results.json').read_text())
            self.assertEqual(failed['status'], 'failed')
            self.assertEqual(list(failed['runs']), ['checked-off'])

    @unittest.skipUnless(shutil.which('zstd'), 'zstd required')
    def test_exact_gate_rejects_a_changed_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / 'a', Path(tmp) / 'b'
            fixture(a)
            shutil.copytree(a, b)
            for root in (a, b):
                shutil.copy2(root / 'simple_checkpoint_0003.bin.zst', root / 'simple_checkpoint_0004.bin.zst')
            target = b / 'simple_checkpoint_0004.bin.zst'
            raw = bytearray(subprocess.check_output(['zstd', '-dc', str(target)]))
            struct.pack_into('<d', raw, 24 + 20 * 160 + 24 + 3 * 41 * 8 + 11 * 80 * 8, 1e-15)
            subprocess.run(['zstd', '-q', '-f', '-o', str(target)], input=raw, check=True)
            self.assertFalse(campaign.exact_gate(a, b, 'simple', [4])['passed'])


if __name__ == '__main__':
    unittest.main()
