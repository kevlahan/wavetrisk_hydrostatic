import copy
import unittest
from compare_complete_gap_195 import check_records,ledger


def record(parts,restart=False,remap=False):
    values=[0.]*67
    for i,v in parts.items():values[i-1]=v
    return {'wall':sum(values),'self':values,'restart':restart,'remap':remap}


class CompleteGapTests(unittest.TestCase):
    def setUp(self):
        self.b={1:{0:record({5:2.,28:1.}),1:record({5:1.,28:1.})},
                2:{0:record({5:3.},remap=True),1:record({5:2.},remap=True)}}
        self.l={1:{0:record({28:1.}),1:record({28:2.})},
                2:{0:record({28:1.},remap=True),1:record({28:2.},remap=True)}}

    def test_signed_partition_conserves(self):
        check_records(self.b,2);check_records(self.l,2)
        result=ledger(self.b,self.l)
        r=result['summaries']['all_ordinary']['rank_mean']
        self.assertAlmostEqual(sum(r['group_differences_seconds'].values()),r['gap_seconds'])
        self.assertLess(r['group_differences_seconds']['Domain boundary handling'],0)
        self.assertEqual(result['summaries']['with_remap']['steps'],1)
        self.assertEqual(result['summaries']['without_remap']['steps'],1)

    def test_actual_slowest_rank_not_component_maximum(self):
        self.b[1][1]=record({5:0.,28:2.5})
        r=ledger(self.b,self.l)['steps'][0]['slowest_rank']
        self.assertEqual(r['selected_ranks'],{'block':0,'legacy':1})
        self.assertEqual(r['wall']['block'],3.)
        self.assertEqual(r['groups']['Domain boundary handling']['block'],1.)

    def test_missing_rank_and_nonconserving_time_fail(self):
        del self.b[1][1]
        with self.assertRaisesRegex(ValueError,'rank coverage'):check_records(self.b,2)
        self.l[1][0]['wall']+=1
        with self.assertRaisesRegex(ValueError,'conserve'):check_records(self.l,2)

    def test_mismatched_sequence_and_events_fail(self):
        other=copy.deepcopy(self.l);del other[2]
        with self.assertRaisesRegex(ValueError,'sequence'):ledger(self.b,other)
        self.l[1][0]['remap']=True
        with self.assertRaisesRegex(ValueError,'flags'):check_records(self.l,2)
        with self.assertRaisesRegex(ValueError,'workload'):ledger(self.b,self.l)

    def test_nonfinite_and_negative_fail(self):
        for value in (float('nan'),float('inf'),-1.):
            other=copy.deepcopy(self.b);other[1][0]['self'][0]=value
            with self.assertRaisesRegex(ValueError,'Invalid time'):check_records(other,2)

    def test_restart_excluded_from_ordinary(self):
        for data in (self.b,self.l):
            for r in data[2].values():r['restart']=True
        out=ledger(self.b,self.l)
        self.assertEqual(out['summaries']['all_ordinary']['steps'],1)
        self.assertNotIn('with_remap',out['summaries'])


if __name__=='__main__':unittest.main()
