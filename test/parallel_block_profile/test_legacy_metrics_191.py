import tempfile
import unittest
from pathlib import Path
from legacy_metrics_191 import collect


def fixture(root):
    for rank in range(2):
        lines=[]
        for i in range(1,66):
            lines += [f'name {i} region {i}',f'region {i} 0 0 0 0 0']
        lines += [f'step 12 0 0 {rank+1} 0']
        lines += [f'self 1 {i} {rank+1 if i==1 else 0} 0' for i in range(1,66)]
        (root/f'legacy-detail-rank-{rank}.txt').write_text('\n'.join(lines)+'\n')


class LegacyMetrics191(unittest.TestCase):
    def test_complete_and_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture(root)
            result=collect(root,2)
            self.assertEqual(result['steps'][0]['wall_mean'],1.5)
            self.assertEqual(result['steps'][0]['max_rank'],1)
            path=root/'legacy-detail-rank-1.txt';original=path.read_text()
            for change in (original.replace('self 1 65 0 0\n',''),
                           original.replace('self 1 1 2 0','self 1 1 0 0'),
                           original.replace('self 1 1 2 0','self 1 1 nan 0'),
                           original+'step 12 0 0 2 0\n'):
                path.write_text(change)
                with self.assertRaises(ValueError):collect(root,2)
            path.unlink()
            with self.assertRaises(ValueError):collect(root,2)

if __name__=='__main__':unittest.main()
