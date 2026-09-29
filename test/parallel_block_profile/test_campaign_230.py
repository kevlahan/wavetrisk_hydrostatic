import unittest
from cluster_230 import campaign
from protocol_230 import screen,confirm
from validation_230 import storage_evidence
from placement_230 import plan


class Campaign(unittest.TestCase):
    def run_case(self,values):
        calls=[];decisions=[]
        def execute(level,n,kind,oracle,phase,suffix):
            calls.append((n,kind,oracle,phase,suffix))
            return values.get((n,kind,phase),1.)
        result=campaign(execute,lambda d:decisions.append(dict(d)))
        return calls,result
    def test_neutral_cleanup_gets_confirmation(self):
        calls,result=self.run_case({})
        self.assertEqual(len(calls),15)
        self.assertEqual(result['confirmation']['decision'],'within_regression_guard')
        self.assertEqual([c[1] for c in calls[-4:]],['baseline','candidate','candidate','baseline'])
        self.assertEqual(sum(c[0]==167 for c in calls),4)
        self.assertFalse(any(c[2] or c[1]=='checked' for c in calls if c[3]!='validation'))
    def test_clear_regression_stops_timing(self):
        calls,result=self.run_case({(16,'candidate','screen'):1.05})
        self.assertEqual(len(calls),9)
        self.assertNotIn('transfer',result)
        self.assertEqual(result['screen']['decision'],'review_regression')
    def test_near_threshold_repetition_bounded(self):
        calls,result=self.run_case({(83,'candidate','confirm'):1.03})
        self.assertEqual(len(calls),19)
        self.assertEqual(result['confirmation']['decision'],'inconclusive')
    def test_disagreement_needs_more_evidence(self):
        batch=list(zip(('baseline','candidate','candidate','baseline'),(1,1.05,1,1)))
        self.assertEqual(confirm([batch])['decision'],'repeat_abba')
    def test_storage_is_collective(self):
        for oracle in (False,True):
            count=167 if oracle else 0
            header='legacy/oracle' if oracle else 'absent (resident block dynamics)'
            text=(f'Domain tendency storage ranks = {count} / 167\nDomain tendency storage: {header}\n')*2
            self.assertEqual(storage_evidence(text,167,'candidate',oracle)['storage_ranks'],[(count,167)]*2)
            with self.assertRaises(ValueError):storage_evidence(text,167,'candidate',not oracle)
            with self.assertRaises(ValueError):storage_evidence(text.replace(f'= {count}', '= 1'),167,'candidate',oracle)
        with self.assertRaises(ValueError):storage_evidence('',167,'checked',False)
    def test_physical_placement(self):
        top=dict(host='test',rows=[dict(cpu=i,socket=i//128,core=i%128,node=i//32,llc=[i//16]) for i in range(256)])
        cases=plan(top)
        self.assertEqual(set(cases),{16,83,167})
        self.assertEqual(cases[167]['socket_rank_counts'],[84,83])
        for n,case in cases.items():self.assertEqual(len({r['cpu'] for r in case['mapping']}),n)

if __name__=='__main__':unittest.main()
