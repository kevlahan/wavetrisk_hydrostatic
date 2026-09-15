import copy
from pathlib import Path
import tempfile
import unittest
from analyze_stage195_campaign import matched_accounting, transactions


class CampaignAccountingTest(unittest.TestCase):
    def runs(self):
        runs={}
        for n in range(1,7):
            for side,wall,phases,slow in [('block',10.,[0,7,1,2,0,0],[0,8,1,3,0,0]),
                                         ('legacy',5.,[0,3,1,1,0,0],[0,4,1,1,0,0])]:
                s={'sequence':1,'restart':0,'remap':0,'wall_mean':wall,'phase_mean':phases,
                   'wall_max':sum(slow),'phase_on_slow_rank':slow}
                checkpoint=copy.deepcopy(s);checkpoint.update(sequence=2,restart=1)
                runs[f'{side}-coarse-{n}']={'workflow':{'steps':[s,checkpoint]}}
        return runs

    def test_additive_slowest_rank_and_restart_exclusion(self):
        a=matched_accounting(self.runs(),'coarse')
        self.assertEqual(a['pairs'][0]['steps'],1)
        self.assertEqual(a['rank_mean']['gap'],5.)
        self.assertEqual(a['slowest_rank']['gap'],6.)
        self.assertEqual(a['slowest_rank']['phase_excess']['adapt_reconstruct'],2.)
        self.assertEqual(a['slowest_rank']['max_abs_residual'],0.)

    def test_unmatched_events_rejected(self):
        r=self.runs();r['block-coarse-1']['workflow']['steps'][0]['remap']=1
        with self.assertRaises(AssertionError):matched_accounting(r,'coarse')

    def test_bad_conservation_rejected(self):
        r=self.runs();r['block-coarse-1']['workflow']['steps'][0]['wall_mean']=12.
        with self.assertRaises(AssertionError):matched_accounting(r,'coarse')

    def test_transaction_overlap_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'workflow195-rank-0.txt'
            p.write_text('step 1 0 0 10\ntransaction 1 1 1 1 4\nevent 1 1 1 1 0 4\nevent 1 2 1 2 3 6\n')
            with self.assertRaises(AssertionError):transactions(Path(d))
            p.write_text('step 1 0 0 10\ntransaction 1 1 1 1 4\nevent 1 1 1 1 0 4\nevent 1 2 1 2 4 6\n')
            self.assertEqual(transactions(Path(d))['seconds']['advance:tag1'],4.)


if __name__=='__main__':unittest.main()
