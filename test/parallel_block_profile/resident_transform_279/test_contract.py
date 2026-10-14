"""Failure-path tests for campaign source, rank coverage and sealing gates."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from common import source_scope,preflight,runtime,J5_HASHES,resource_summary,summarize,ORDER
from contract import evidence
from campaign_gate import check,BASE

class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def rows(self):
        for rank in range(2):
            lines=[]
            for step in range(1,49):
                values=[step,step,1,1,1-rank,1,100*step,80*step,60*step,step,20*step,40*step,10*step,step,30*step,20*step,step*.1,step*.01]
                lines.append(' '.join(map(str,values)))
            (self.root/f'rhs278-cost-rank-{rank}.txt').write_text('1 1 10 20 30 30 12 12 1 8 8 1 100 200 300 400 400 0.01\n')
            (self.root/f'rhs278-production-rank-{rank}.txt').write_text('\n'.join(lines)+'\n')
    def test_above_domain_rank_does_real_work(self):
        self.rows();result=evidence(self.root,2,True,domains=1)
        self.assertEqual(result['zero_domain_ranks'],1)
        self.assertEqual(result['zero_domain_rank_primitive_events'],4800)
    def test_reject_missing_rank(self):
        self.rows();(self.root/'rhs278-production-rank-1.txt').unlink()
        with self.assertRaises(ValueError):evidence(self.root,2,True,domains=1)
    def test_reject_missing_rk_stage(self):
        self.rows();p=self.root/'rhs278-production-rank-1.txt';p.write_text('\n'.join(p.read_text().splitlines()[:-1]))
        with self.assertRaises(ValueError):evidence(self.root,2,True,domains=1)
    def test_reject_idle_zero_domain_owner(self):
        self.rows();p=self.root/'rhs278-production-rank-1.txt';rows=[x.split() for x in p.read_text().splitlines()]
        for row in rows:row[6]='0'
        p.write_text('\n'.join(' '.join(x) for x in rows))
        with self.assertRaises(ValueError):evidence(self.root,2,True,domains=1)
    def test_reject_cost_epoch_gap(self):
        self.rows();p=self.root/'rhs278-cost-rank-1.txt';p.write_text(p.read_text().replace('1 1 10','1 2 10'))
        with self.assertRaisesRegex(ValueError,'cost epoch coverage'):evidence(self.root,2,True,domains=1)
    def test_bare_control_rejects_audit(self):
        self.rows()
        with self.assertRaises(ValueError):evidence(self.root,2,False)

if __name__=='__main__':unittest.main()
