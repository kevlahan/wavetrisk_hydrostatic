"""Fail closed on incomplete evidence; keep a failed residency hypothesis visible."""
import tempfile,unittest
from pathlib import Path
from contract import evidence

class Evidence(unittest.TestCase):
    def make(self,root,diff=0,ref=0):
        p=root/'primitive276-rank-0.txt'
        counts=[1,1,4,1,1,1,0,14,1,1,diff,0,0,9,diff,ref]
        lines=['RANK 0 1 40 40']
        for transaction in (1,2,3,9,17,33,41,48):
            for k in (1,30):
                counts[0:2]=[transaction,k]
                lines.append('PHASE '+' '.join(map(str,counts+[.01]*7)))
        p.write_text('\n'.join(lines)+'\n')
        return p
    def test_valid_and_failed_residency_remain_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.make(root)
            result=evidence(root,1,True,40)
            self.assertTrue(result['resident_gate_passed']);self.assertFalse(result['performance_qualified'])
            self.make(root,diff=1)
            self.assertFalse(evidence(root,1,True,40)['resident_gate_passed'])
    def test_bad_reference_cannot_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.make(root,ref=1)
            with self.assertRaises(ValueError):evidence(root,1,True,40)
    def test_missing_duplicate_and_truncated_evidence_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=self.make(root)
            with self.assertRaises(ValueError):evidence(root,2,True,40)
            p.write_text(p.read_text()+p.read_text().splitlines()[1]+'\n')
            with self.assertRaises(ValueError):evidence(root,1,True,40)
            p=self.make(root);p.write_text(p.read_text().replace(' 9 ',' 8 '))
            with self.assertRaises(ValueError):evidence(root,1,True,40)
    def test_control_must_not_contain_prototype_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            self.assertEqual(evidence(root,1,False),dict(enabled=False))
            self.make(root)
            with self.assertRaises(ValueError):evidence(root,1,False)

if __name__=='__main__':unittest.main()
