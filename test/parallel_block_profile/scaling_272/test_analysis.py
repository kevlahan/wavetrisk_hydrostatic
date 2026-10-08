"""Exercise attribution using synthetic rank identities, zero owners and corrupted payload counts."""
import copy
import unittest
from cost_analysis import ownership_groups, INCLUDED
from phases import NAMES, IDS


def records():
    ranks=[]
    for r in range(2):
        steps=[]
        for s in range(12):
            rows=[dict(calls=0,self_wall=0.,inclusive_wall=0.,work=[0,0,0]) for _ in NAMES]
            rows[IDS['OWNERS']-1].update(calls=4,work=[4 if r==1 else 0,8,4])
            rows[IDS['PATCHES']-1]['work']=[40 if r==1 else 0,20,4]
            rows[IDS['BASIC']-1].update(self_wall=.2 if r==1 else .001,work=[80 if r==1 else 0,0,0])
            rows[IDS['CAPTURE_MPI']-1].update(self_wall=.001 if r==1 else .2,work=[4,12,12])
            steps.append(rows)
        ranks.append(steps)
    return ranks


class Analysis(unittest.TestCase):
    def test_ownership_uses_census_and_preserves_distributions(self):
        result=ownership_groups(records(),1)
        self.assertEqual(result['Domain_owners']['ranks_by_step']['2'],[1])
        self.assertEqual(result['without_Domains']['ranks_by_step']['2'],[0])
        self.assertEqual(result['Domain_owners']['rank_step_samples'],9)
        self.assertEqual(result['Domain_owners']['phases']['Domain basic operators']['self_wall']['median'],.2)
        self.assertEqual(result['without_Domains']['phases']['scalar capture MPI']['self_wall']['maximum'],.2)

    def test_incomplete_or_inconsistent_census_rejected(self):
        for change in ('missing','global','impossible_work'):
            rows=records()
            if change=='missing':rows[0][1][IDS['OWNERS']-1]['work'][2]=3
            if change=='global':rows[1][1][IDS['OWNERS']-1]['work'][0]=8
            if change=='impossible_work':rows[0][1][IDS['BASIC']-1]['work'][0]=1
            with self.assertRaises(ValueError):ownership_groups(rows,1)

    def test_remote_payload_corruption_rejected(self):
        rows=records();rows[1][1][IDS['CAPTURE_MPI']-1]['work'][1]+=1
        with self.assertRaises(ValueError):ownership_groups(rows,1)

if __name__=='__main__':unittest.main()
