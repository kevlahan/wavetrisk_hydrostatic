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
import cluster_196 as campaign
from test_cluster_187 import fixture, LOG


class Cluster196(unittest.TestCase):
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
                    'print("Stage 193 requested scalar inverse columns = 1")\n'
                    'print("Stage 194 requested inner inverse geometry = "+os.environ["WAVETRISK_STAGE194_INNER_GEOMETRY"])\n')
                binary.chmod(0o755);binaries[kind]=binary
            args=argparse.Namespace(**binaries,fixture=seed,out=root/'results',ranks=83,repeats=1,profiles=True,
                                    expected_checkpoints=[4],require_post_restart_remap=True,build_notes='test')
            env={'PATH':str(root)+os.pathsep+os.environ['PATH'],'SLURM_JOB_ID':'test','SLURM_NTASKS':'83'}
            for failure in ('none','workload','launch'):
                fail=failure!='none'
                args.out=root/failure
                with patch.dict(os.environ,{**env,'INJECT_WORKLOAD':str(int(failure=='workload')),
                                             'INJECT_EXIT':str(int(failure=='launch'))}), \
                     patch('cluster_196.snapshot_allocation',return_value={'ranks':83}), \
                     patch('cluster_196.node_memory',return_value=[]), \
                     patch('cluster_196.collect_legacy',return_value={'ranks':83}), \
                     patch('cluster_196.stage191_evidence',return_value={}), patch('cluster_196.collect_workflow',side_effect=lambda root,ranks,mode: dict(mode=mode,ordinary_wall_max_median=1.,ordinary_rank_mean_wall=1.,ordinary_phase_mean={n:1/6 for n in ('setup','advance','physics','adapt_reconstruct','remap','output_restart')},steps=[dict(sequence=1,restart=0,remap=0,transaction_mean={'1':{'calls_global':5},'2':{'calls_global':4}}),dict(sequence=2,restart=1,remap=1,transaction_mean={'1':{'calls_global':5},'2':{'calls_global':4}})])),patch('cluster_196.collect_operations',return_value={'ordinary_context_wall_means':{'WT':0.5,'outside_WT':0.5},'ordinary_context_group_means':{c:{g:0. for g in __import__('analyze_operations_196').GROUPS} for c in ('WT','outside_WT')}}),contextlib.redirect_stdout(io.StringIO()):
                    if failure=='launch':
                        with self.assertRaisesRegex(RuntimeError,'checked-off exited 74'):campaign.run(args)
                    elif fail:
                        with self.assertRaisesRegex(ValueError,'Workload differs'):campaign.run(args)
                    else:campaign.run(args)
                result=json.loads((args.out/'results.json').read_text())
                self.assertEqual(result['status'],'failed' if fail else 'passed')
                if not fail:
                    self.assertEqual(len(result['runs']),18)
                    self.assertEqual(len(result['gates']),46)
                    self.assertEqual(result['runs']['block-minimal-1']['stage194_evidence']['requested_inner_geometry'],'1')
                    self.assertEqual(result['runs']['block-coarse-1']['environment']['WAVETRISK_STAGE191_MASS_OVERLAP'],'0')
                    self.assertTrue(result['gates']['cross-mode-candidate-warmup']['passed'])
                    self.assertTrue(result['gates']['legacy-screen-legacy-bare-1']['passed'])
                    self.assertEqual(result['summary']['modes']['bare']['legacy']['repeats'],1)
                else:self.assertNotIn('baseline-1',result['runs'])

if __name__=='__main__':unittest.main()
