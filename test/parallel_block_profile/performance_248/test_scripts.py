"""Offline checks of campaign selection, physical-core maps and archive integrity."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

PERF=Path(__file__).resolve().parent
sys.path.insert(0,str(PERF))
sys.path.insert(0,str(PERF.parents[0]))
from build import commands
from collect import collect, verify
from placement import plan
from timing import screen, confirm


class ScriptChecks(unittest.TestCase):
    def test_matched_builds(self):
        kit=Path(__file__).resolve().parents[3]/'.stage247-validation/bbserv-candidate2'
        for level in ('J4','J5'):
            for comparison in ('mirrored','feature-off'):
                configs=commands(kit,level,comparison)
                for key in ('candidate','baseline'):
                    self.assertIn('PARAM=param_'+level,configs[key])
                    self.assertEqual(configs[key][:2],['make','-j1'])
                self.assertIn('MIRROR_FREE=on',configs['candidate'])
                self.assertIn('MIRROR_FREE=off',configs['baseline'])
                self.assertIn('BLOCK_RK='+('off' if comparison=='feature-off' else 'on'),configs['baseline'])
                self.assertIn('COEFFICIENT_NATIVE='+('off' if comparison=='feature-off' else 'on'),configs['baseline'])

    def test_physical_core_maps(self):
        rows=[dict(cpu=s*128+i,core=i,socket=s,node=s,llc=list(range(s*128,(s+1)*128)))
              for s in (0,1) for i in range(128)]
        top=dict(rows=rows,host='bb13')
        for n in (16,41,83,167):
            case=plan(top,n)
            self.assertEqual(len(case['mapping']),n)
            self.assertEqual(len({r['cpu'] for r in case['mapping']}),n)
            self.assertEqual(case['socket_rank_counts'],[n,0] if n==16 else [(n+1)//2,n//2])

    def test_run_level_decision(self):
        self.assertEqual(screen(1,1.04)['decision'],'review_regression')
        batch=[('baseline',1),('candidate',.95),('candidate',.96),('baseline',1.01)]
        result=confirm([batch])
        self.assertEqual(result['decision'],'within_regression_guard')
        self.assertEqual(result['independent_runs_per_version'],2)

    def test_review_archive(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)/'results';root.mkdir()
            (root/'results.json').write_text(json.dumps(dict(status='complete')))
            (root/'run.log').write_text('sample\n')
            archive=Path(temp)/'review.tar.gz'
            collect(root,archive)
            self.assertEqual(verify(archive)['verified_files'],2)
            archive.write_bytes(archive.read_bytes()[:20])
            with self.assertRaises(Exception):verify(archive)


if __name__=='__main__':unittest.main()
