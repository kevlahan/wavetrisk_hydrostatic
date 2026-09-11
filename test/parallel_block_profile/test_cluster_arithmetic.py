"""Cross-mode gates must reject off/check differences, not default/off differences."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import build_arithmetic_188 as builder
import cluster_187 as campaign
from test_cluster_187 import fixture, LOG


class ClusterArithmeticTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('zstd'), 'zstd required')
    def test_default_difference_is_diagnostic_but_crossmode_difference_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root / 'fixture'
            fixture(seed)
            launch = root / 'srun'
            launch.write_text('#!/usr/bin/env python3\nimport os,sys\n'
                              'i=sys.argv.index("./climateJ5");os.execv(sys.argv[i],sys.argv[i:])\n')
            launch.chmod(0o755)
            binaries = {}
            def executable(kind, value):
                path = root / kind
                path.write_text('#!/usr/bin/env python3\nimport subprocess,struct\n'
                    'raw=bytearray(subprocess.check_output(["zstd","-dc","simple_checkpoint_0003.bin.zst"]))\n'
                    f'struct.pack_into("<d",raw,24+20*160+24+3*41*8+11*80*8,{value!r})\n'
                    'subprocess.run(["zstd","-q","-o","simple_checkpoint_0004.bin.zst"],input=raw,check=True)\n'
                    'print(' + repr(LOG.replace('Restarting from checkpoint 4',
                        'Restarting from checkpoint 4\nRemapping vertical coordinates')) + ')\n')
                path.chmod(0o755)
                return path
            for kind in ('baseline', 'candidate', 'checked'):
                binaries[kind] = executable(kind, 1.0 if kind == 'baseline' else 0.0)
            args = argparse.Namespace(**binaries, fixture=seed, out=root / 'passed', ranks=83,
                                      repeats=1, profiles=False, expected_checkpoints=[4],
                                      require_post_restart_remap=True, build_notes='Synthetic arithmetic control',
                                      arithmetic_control=True)
            env = {'PATH': str(root) + os.pathsep + os.environ['PATH'],
                   'SLURM_JOB_ID': 'test', 'SLURM_NTASKS': '83'}
            with patch.dict(os.environ, env), contextlib.redirect_stdout(io.StringIO()):
                campaign.run(args)
            result = json.loads((args.out / 'results.json').read_text())
            self.assertEqual(result['status'], 'passed')
            self.assertFalse(result['default_vs_off_diagnostic']['passed'])
            self.assertTrue(all(g['passed'] for g in result['gates'].values()))
            self.assertTrue(result['gates']['candidate-1']['reference'].endswith('checked-off'))
            self.assertTrue(result['gates']['baseline-1']['reference'].endswith('baseline-warmup'))
            executable('candidate', 1e-15)
            args.out = root / 'failed'
            with patch.dict(os.environ, env), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, 'Exact checkpoint gate failed for candidate-warmup'):
                    campaign.run(args)
            failed = json.loads((args.out / 'results.json').read_text())
            self.assertEqual(failed['status'], 'failed')
            self.assertNotIn('baseline-1', failed['runs'])

    def test_build_rejects_changed_source_and_existing_build_before_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, tree = root / 'kit', root / 'tree'
            kit.mkdir(); tree.mkdir()
            (tree / 'source.f90').write_text('source')
            manifest = {'sources': {'source.f90': 'wrong'}}
            target = kit / 'arithmetic-188-sources.json'
            target.write_text(json.dumps(manifest))
            with patch.object(builder.subprocess, 'run') as invoke:
                with self.assertRaisesRegex(ValueError, 'Source differs'):
                    builder.build(kit, tree)
                manifest['sources']['source.f90'] = builder.digest(tree / 'source.f90')
                target.write_text(json.dumps(manifest))
                (tree / 'build-188-default').mkdir()
                with self.assertRaisesRegex(ValueError, 'Fresh build directory required'):
                    builder.build(kit, tree)
                invoke.assert_not_called()


if __name__ == '__main__':
    unittest.main()
