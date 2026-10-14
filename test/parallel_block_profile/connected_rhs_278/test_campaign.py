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
    def test_fail_closed_without_active_science(self):
        gate=check(self.root,{'accepted_source':BASE})
        self.assertFalse(gate['sealing_ready']);self.assertFalse(gate['launch_authorized'])
    def test_unsealed_kit_cannot_launch(self):
        with patch('common.verify'),patch('common.verify_packaged_grid'):
            with self.assertRaisesRegex(ValueError,'Unsealed campaign'):
                preflight(self.root,self.root)
    def test_source_overlay_rejects_extra_change(self):
        a={'Makefile':'a','src/ops.f90':'b'};b={'Makefile':'c','src/ops.f90':'b'}
        (self.root/'source-files.json').write_text(json.dumps(dict(stage266=a,main=a,candidate=b)))
        (self.root/'overlay-files.json').write_text(json.dumps(['Makefile']))
        identity=dict(changed_source=['Makefile'],candidate_source_sha256={'Makefile':'c'},baseline_source_sha256={'Makefile':'a'})
        (self.root/'IDENTITY.json').write_text(json.dumps(identity))
        self.assertEqual(source_scope(self.root),['Makefile'])
        b['src/ops.f90']='unexpected'
        (self.root/'source-files.json').write_text(json.dumps(dict(stage266=a,main=a,candidate=b)))
        with self.assertRaises(ValueError):source_scope(self.root)
    def test_resource_summary_requires_actual_child_measurements(self):
        row=dict(returncode=0,resource_collector='os.wait4',max_rss_kib=100,major_faults=0,minor_faults=2,user_cpu_seconds=1.,system_cpu_seconds=.5)
        self.assertEqual(resource_summary([row])['solver_cpu_seconds'],1.5)
        row['max_rss_kib']=0
        with self.assertRaises(ValueError):resource_summary([row])
    def test_dual_goal_report_keeps_gap_and_scaling_separate(self):
        rows=[dict(ranks=n,kind=k,ordinary={'median':(1. if k=='main' else 2.)},trajectory_step_seconds=24.,resources={}) for n,k in ORDER]
        result=summarize(rows)
        self.assertEqual(result['by_ranks']['160']['candidate_minus_main_percent'],100.)
        self.assertEqual(result['by_ranks']['160']['reduction_needed_for_main_parity_percent'],50.)
        self.assertEqual(result['speedup_160_to_256']['candidate'],1.)
        self.assertEqual(result['relative_efficiency_160_to_256_percent']['candidate'],62.5)
        self.assertFalse(result['automatic_promotion'])
    def test_grid_identity_retained(self):
        self.assertEqual(J5_HASHES['grids/HR95JT_WT_004'],'d4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d')

if __name__=='__main__':unittest.main()
