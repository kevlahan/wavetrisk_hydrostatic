#!/usr/bin/env python3
"""Numerically meaningful scheduler/statistics tests; no cluster required."""
import unittest
from pathlib import Path
from protocol_215 import ordinary,screen,confirm
from cluster_215 import campaign,cases_for
from analyze import analyze
import tempfile

class Protocol(unittest.TestCase):
    def test_startup_checkpoint_reload_exclusions(self):
        lines=['Restarting from checkpoint 3']
        for i in range(10):
            if i==4:lines.append('Saving checkpoint')
            if i==5:lines.append('Restarting from checkpoint 4')
            if i==7:lines.append('Remapping vertical coordinates')
            lines.append(f'{i+1} 0.3 d dt = 1 cpu = {100 if i in (0,4,5) else 1+i/100}')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run.log';p.write_text('\n'.join(lines));a=analyze(p);r=ordinary(p.read_text(),a)
        self.assertEqual([r['index'] for r in r['excluded']],[0,4,5])
        self.assertEqual(len(r['included']),7);self.assertAlmostEqual(r['median'],1.06)
        self.assertIn(7,[v['index'] for v in r['included']])

    def test_reload_inside_checkpoint_transaction(self):
        lines=['Restarting from checkpoint 3']
        for i in range(12):
            if i==7:lines += ['Saving checkpoint 4','Restarting from checkpoint 4','Saving field 4']
            if i==9:lines.append('Remapping vertical coordinates')
            lines.append(f'{i+1} 0.3 d dt = 1 cpu = {100 if i in (0,7,8) else 1}')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run.log';p.write_text('\n'.join(lines));r=ordinary(p.read_text(),analyze(p))
        self.assertEqual([v['index'] for v in r['excluded']],[0,7,8])
        self.assertEqual(len(r['included']),9);self.assertEqual(r['median'],1)
        self.assertIn(9,[v['index'] for v in r['included']])

    def test_checkpoint_without_reload_does_not_exclude_next_step(self):
        lines=[]
        for i in range(10):
            if i==4:lines.append('Saving checkpoint 4')
            lines.append(f'{i+1} 0.3 d dt = 1 cpu = 1')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run.log';p.write_text('\n'.join(lines));r=ordinary(p.read_text(),analyze(p))
        self.assertEqual([v['index'] for v in r['excluded']],[0,4])

    def run_schedule(self,small=(2,1.8),large=(1,.8),confirmation=None):
        seen=[];series=iter(confirmation or [1,.8,.8,1]);decisions=[]
        def execute(case,kind,phase,suffix,oracle):
            seen.append((case,kind,phase,suffix,oracle))
            if phase=='validation':return 999 # Must not enter timing estimates.
            if phase in ('screen','transfer'):return (small if phase=='screen' else large)[kind=='candidate']
            return next(series)
        result=campaign(execute,lambda d:decisions.append(dict(d)))
        return result,seen

    def test_16_regression_stops(self):
        r,s=self.run_schedule(small=(1,1.1));self.assertEqual(len(s),6)
        self.assertEqual(r['screen']['decision'],'reject_regression')
        self.assertEqual({row[0] for row in s},{'n16-s1'})

    def test_flat_stops_without_repeating(self):
        r,s=self.run_schedule(small=(1,.99));self.assertEqual(len(s),6)
        self.assertEqual(r['screen']['decision'],'stop_no_clear_gain')

    def test_83_transfer_precedes_confirmation(self):
        r,s=self.run_schedule(large=(1,1.2));self.assertEqual(len(s),8)
        self.assertNotIn('confirmation',r)

    def test_small_83_gain_does_not_launch_confirmation(self):
        r,s=self.run_schedule(large=(1,.95));self.assertEqual(len(s),8)
        self.assertEqual(r['transfer']['decision'],'stop_below_target')

    def test_data_save_excluded_independently(self):
        log='\n'.join(('Saving field' if i==3 else '')+'\n'+f'{i} 0.3 d dt = 1 cpu = 1' for i in range(9))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run.log';p.write_text(log);r=ordinary(log,analyze(p))
        self.assertEqual([v['index'] for v in r['excluded']],[0,3])

    def test_fresh_abba(self):
        r,s=self.run_schedule();self.assertEqual(len(s),12)
        self.assertEqual([v[1] for v in s[-4:]],['baseline','candidate','candidate','baseline'])
        self.assertEqual(r['confirmation']['independent_runs_per_version'],2)
        self.assertEqual(r['confirmation']['decision'],'meets_target')
        self.assertAlmostEqual(r['confirmation']['gain_seconds'],.2)

    def test_only_near_or_disagreement_repeat(self):
        r,s=self.run_schedule(confirmation=[1,.9,.9,1]*2)
        self.assertEqual(len(s),16);self.assertEqual(r['confirmation']['decision'],'inconclusive')
        r,s=self.run_schedule(confirmation=[1,.95,.85,1]*2)
        self.assertEqual(len(s),16);self.assertTrue(r['confirmation']['disagreement'])
        r,s=self.run_schedule(confirmation=[1,.95,.95,1]);self.assertEqual(len(s),12)
        self.assertEqual(r['confirmation']['decision'],'below_target')

    def test_validation_failure_stops(self):
        seen=[]
        def fail(*a):seen.append(a);raise ValueError('Numerical gate')
        with self.assertRaises(ValueError):campaign(fail,lambda _:None)
        self.assertEqual(len(seen),1)

    def test_topology_is_fixed_physical_core_map(self):
        rows=[dict(cpu=s*128+i,core=s*128+i,socket=s,node=s,llc=list(range(s*128+(i//8)*8,s*128+(i//8+1)*8))) for s in range(2) for i in range(128)]
        cases=cases_for(dict(rows=rows,host='synthetic'))
        self.assertEqual(set(cases),{'n16-s1','n83-s2'})
        self.assertEqual(len(cases['n16-s1']['mapping']),16)
        self.assertEqual(len({tuple(r['llc']) for r in cases['n16-s1']['mapping']}),16)
        self.assertEqual([sum(r['socket']==s for r in cases['n83-s2']['mapping']) for s in range(2)],[42,41])

    def test_bad_confirmation_order_rejected(self):
        with self.assertRaises(ValueError):confirm([[('baseline',1),('candidate',.8),('baseline',1),('candidate',.8)]])

if __name__=='__main__':unittest.main()
