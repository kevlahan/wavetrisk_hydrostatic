import tempfile,unittest
from pathlib import Path
from prebound_contract import evidence
class CompactPublication(unittest.TestCase):
    def fixture(self,p):
        (p/'run.log').write_text('Transform280 publication: 0 10 2 6\n'*12)
        (p/'transform279-cost-rank-0.txt').write_text('1 1 1 1 0 0 10 1 1 0 0 1 1 0.5\n'*12)
    def test_current_bindings_and_timing_isolation(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);r=evidence(p,1,True);self.assertEqual(r['avoided_numerical_rows'],96)
            with self.assertRaises(ValueError):evidence(p,1,False)
    def test_wrong_counts_missing_epoch_and_impossible_extent_rejected(self):
        for mode in ('census','missing','extent'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as t:
                p=Path(t);self.fixture(p)
                if mode=='census':(p/'run.log').write_text('Transform280 publication: 0 11 2 6\n'*12)
                elif mode=='extent':(p/'run.log').write_text('Transform280 publication: 0 10 5 6\n'*12)
                else:(p/'run.log').write_text('Transform280 publication: 0 10 2 6\n'*11)
                with self.assertRaises(ValueError):evidence(p,1,True)
