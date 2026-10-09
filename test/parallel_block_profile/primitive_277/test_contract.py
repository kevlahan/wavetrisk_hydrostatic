import tempfile,unittest
from pathlib import Path
from contract import evidence

class Contract(unittest.TestCase):
    def fixture(self,root):
        (root/'run.log').write_text('Stage277 cached fused primitive = on; assignment checks = T\n')
        for r in range(4):
            rows=[f'RANK {r} 4 3 {int(r<3)} 30 4']
            for t in range(1,49):
                n=64 if r<3 else 0
                rows.append(' '.join(map(str,['RHS',t,1,1,n*30*19]+[n]*30)))
            (root/f'primitive277-rank-{r}.txt').write_text('\n'.join(rows)+'\n')
    def test_full_and_zero_domain_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p);r=evidence(p,4,True,3)
            self.assertEqual(r['nodes'],3*48*30*64)
            self.assertEqual(r['zero_domain_ranks'],1)
            self.assertEqual(r['zero_domain_rank_nodes'],0)
    def test_missing_layer_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p);f=p/'primitive277-rank-0.txt';s=f.read_text().splitlines();row=s[1].split();row[-1]='0';row[4]=str(19*29*64);s[1]=' '.join(row);f.write_text('\n'.join(s))
            with self.assertRaises(ValueError):evidence(p,4,True,3)
    def test_missing_rhs_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p);f=p/'primitive277-rank-0.txt';f.write_text('\n'.join(f.read_text().splitlines()[:-1]))
            with self.assertRaises(ValueError):evidence(p,4,True,3)
    def test_checks_forbidden_in_timing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p)
            with self.assertRaises(ValueError):evidence(p,4,False,3)
            for f in p.glob('primitive277-rank-*'):f.unlink()
            with self.assertRaises(ValueError):evidence(p,4,False,3)
            (p/'run.log').write_text('Stage277 cached fused primitive = on; assignment checks = F\n')
            self.assertFalse(evidence(p,4,False,3)['enabled'])
if __name__=='__main__':unittest.main()
