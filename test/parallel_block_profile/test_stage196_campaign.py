"""Regression checks for archive-versus-raw comparison acceptance."""
import copy
import unittest
from analyze_stage196_campaign import close

class ArchiveComparisonTest(unittest.TestCase):
    def setUp(self):
        self.value={'rank':82,'steps':[{'seconds':0.1234567890123,'calls':177,'valid':True}]}

    def test_float_roundoff_is_allowed(self):
        candidate=copy.deepcopy(self.value)
        candidate['steps'][0]['seconds']+=1e-14
        close(self.value,candidate)

    def test_timing_corruption_is_rejected(self):
        candidate=copy.deepcopy(self.value)
        candidate['steps'][0]['seconds']+=1e-6
        with self.assertRaises(AssertionError):close(self.value,candidate)

    def test_missing_rank_or_step_is_rejected(self):
        for key in ('rank','steps'):
            candidate=copy.deepcopy(self.value);del candidate[key]
            with self.assertRaises(AssertionError):close(self.value,candidate)
        candidate=copy.deepcopy(self.value);candidate['steps']=[]
        with self.assertRaises(AssertionError):close(self.value,candidate)

    def test_counts_and_status_are_exact(self):
        for key,val in (('calls',178),('valid',False)):
            candidate=copy.deepcopy(self.value);candidate['steps'][0][key]=val
            with self.assertRaises(AssertionError):close(self.value,candidate)

    def test_nonfinite_timing_is_rejected(self):
        for val in (float('nan'),float('inf')):
            candidate=copy.deepcopy(self.value);candidate['steps'][0]['seconds']=val
            with self.assertRaises(AssertionError):close(self.value,candidate)

if __name__=='__main__':unittest.main()
