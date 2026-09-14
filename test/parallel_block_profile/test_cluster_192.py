"""Exercise three-way scheduling, numerical/workload gates and fail-closed runs."""
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
import cluster_192 as campaign
from test_cluster_187 import fixture, LOG


class Cluster192(unittest.TestCase):
    def test_missing_or_invalid_build_notes_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'builds.json'
            args = argparse.Namespace(build_notes_file=path)
            with self.assertRaises(FileNotFoundError): campaign.read_build_notes(args)
            path.write_text('invalid json')
            with self.assertRaises(ValueError): campaign.read_build_notes(args)
            path.write_text('[]')
            with self.assertRaises(ValueError): campaign.read_build_notes(args)

    def test_selector_and_production_installation_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'run.log').write_text('Stage 191 requested install/mass-overlap = 1 1\n')
            counters={k+'  ':0 for k in ('source proof payload bytes','receive proof payload bytes',
                'copy proof payload bytes','install copied payload bytes')}
            counters['install moved payload bytes  ']=100
            analysis={'detail_windows':[{'counters':counters}]}
            self.assertEqual(campaign.stage191_evidence(root,analysis,'1','1',True)['install moved payload bytes'],100)
            with self.assertRaises(ValueError):campaign.stage191_evidence(root,analysis,'0','1',False)
            counters['copy proof payload bytes  ']=1
            with self.assertRaises(ValueError):campaign.stage191_evidence(root,analysis,'1','1',True)

    def test_rotation_and_summary(self):
        jobs=campaign.schedule(3,True)
        timed=[j[1] for j in jobs if j[4]=='timing']
        self.assertEqual(timed,['baseline','candidate','legacy','candidate','legacy','baseline','legacy','baseline','candidate'])
        runs={str(i):dict(kind=k,phase='timing',elapsed_seconds={'baseline':4,'candidate':3,'legacy':1}[k],
                         ordinary_step_median={'baseline':2,'candidate':1.5,'legacy':.5}[k]) for i,k in enumerate(timed)}
        runs['warmup']=dict(kind='candidate',phase='warmup',elapsed_seconds=999,ordinary_step_median=999)
        s=campaign.summarize(runs)
        self.assertEqual(s['observed_elapsed_reduction_percent'],25)
        self.assertEqual(s['block_vs_legacy']['candidate']['ordinary_step_ratio'],3)

    @unittest.skipUnless(shutil.which('zstd'),'zstd required')
    def test_legacy_screen_is_not_exact_and_rejects_large_difference(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a',Path(tmp)/'b';fixture(a);shutil.copytree(a,b)
            for root in (a,b):shutil.copy2(root/'simple_checkpoint_0003.bin.zst',root/'simple_checkpoint_0004.bin.zst')
            target=b/'simple_checkpoint_0004.bin.zst'
            raw=bytearray(subprocess.check_output(['zstd','-dc',str(target)]))
            offset=24+20*160+24+3*41*8+11*80*8
            for delta,passed in ((1e-15,True),(1e-5,False)):
                struct.pack_into('<d',raw,offset,delta)
                subprocess.run(['zstd','-q','-f','-o',str(target)],input=raw,check=True)
                self.assertEqual(campaign.legacy_gate(a,b,'simple',[4])['passed'],passed)
                self.assertFalse(campaign.exact_gate(a,b,'simple',[4])['passed'])

    @unittest.skipUnless(shutil.which('zstd'),'zstd required')
    def test_campaign_and_workload_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);seed=root/'fixture';fixture(seed)
            launch=root/'srun'
            launch.write_text('#!/usr/bin/env python3\nimport os,sys\ni=sys.argv.index("./climateJ5")\nos.execv(sys.argv[i],sys.argv[i:])\n')
            launch.chmod(0o755)
            binaries={}
            log=LOG.replace('Restarting from checkpoint 4','Restarting from checkpoint 4\nRemapping vertical coordinates')
            for kind in ('baseline','candidate','checked','legacy','legacy_profile'):
                binary=root/kind
                binary.write_text('#!/usr/bin/env python3\nfrom pathlib import Path\nimport shutil,json,os\n'
                    'if os.environ.get("INJECT_EXIT")=="1": raise SystemExit(74)\n'
                    'shutil.copy2("simple_checkpoint_0003.bin.zst","simple_checkpoint_0004.bin.zst")\n'
                    'for i in range(83):\n'
                    ' Path(f"rank-{i:04d}.json").write_text(json.dumps(dict(rank=i)))\n'
                    ' Path(f"rank-{i:04d}.resources").write_text("Maximum resident set size (kbytes): 100\\nMajor (requiring I/O) page faults: 0\\nMinor (reclaiming a frame) page faults: 1\\n")\n'
                    'text='+repr(log)+'\n'
                    'if os.environ.get("INJECT_WORKLOAD")=="1" and Path.cwd().name=="legacy-warmup": text=text.replace("dof = 999","dof = 998")\n'
                    'print(text)\n'
                    'print("Stage 192 requested sparse outer install = "+os.environ["WAVETRISK_STAGE192_SPARSE_OUTER"])\n')
                binary.chmod(0o755);binaries[kind]=binary
            args=argparse.Namespace(**binaries,fixture=seed,out=root/'results',ranks=83,repeats=1,profiles=True,
                                    expected_checkpoints=[4],require_post_restart_remap=True,build_notes='test')
            env={'PATH':str(root)+os.pathsep+os.environ['PATH'],'SLURM_JOB_ID':'test','SLURM_NTASKS':'83'}
            for failure in ('none','workload','launch'):
                fail=failure!='none'
                args.out=root/failure
                with patch.dict(os.environ,{**env,'INJECT_WORKLOAD':str(int(failure=='workload')),
                                             'INJECT_EXIT':str(int(failure=='launch'))}), \
                     patch('cluster_192.snapshot_allocation',return_value={'ranks':83}), \
                     patch('cluster_192.node_memory',return_value=[]), \
                     patch('cluster_192.collect_legacy',return_value={'ranks':83}), \
                     patch('cluster_192.stage191_evidence',return_value={}), patch('analyze_inverse_192.compare_profiles',return_value={}),contextlib.redirect_stdout(io.StringIO()):
                    if failure=='launch':
                        with self.assertRaisesRegex(RuntimeError,'checked-off exited 74'):campaign.run(args)
                    elif fail:
                        with self.assertRaisesRegex(ValueError,'Workload differs'):campaign.run(args)
                    else:campaign.run(args)
                result=json.loads((args.out/'results.json').read_text())
                self.assertEqual(result['status'],'failed' if fail else 'passed')
                if not fail:
                    self.assertEqual(len(result['runs']),12)
                    self.assertEqual(result['runs']['candidate-control']['stage192_evidence']['requested_sparse_outer'],'0')
                    self.assertEqual(result['runs']['candidate-1']['environment']['WAVETRISK_STAGE191_MASS_OVERLAP'],'0')
                    self.assertTrue(result['gates']['cross-mode-candidate-warmup']['passed'])
                    self.assertTrue(result['gates']['legacy-screen-legacy-1']['passed'])
                    self.assertEqual(len(result['summary']['legacy']['elapsed_seconds']),1)
                else:self.assertNotIn('baseline-1',result['runs'])

if __name__=='__main__':unittest.main()
