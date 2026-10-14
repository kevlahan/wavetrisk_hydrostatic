import unittest,tempfile
from pathlib import Path
from transform_contract import evidence,check_global_census
from diagnostic import environment,region

class Transform(unittest.TestCase):
    def fixture(self,root,oracle=False):
        for rank in range(2):
            row=lambda k:[k,1,1,1-rank,1,3,4,5,int(oracle),100,100,1,27,11,22,1,1,30,40]
            (root/f'transform279-production-rank-{rank}.txt').write_text('\n'.join(' '.join(map(str,row(k))) for k in range(1,37)))
            (root/f'transform279-cost-rank-{rank}.txt').write_text('1 1 3 2 0 0 2 11 22 0 0 1 1 0.5\n')
            (root/f'transform279-work-rank-{rank}.txt').write_text('1 1 1 2 3 4 5\n')
    def test_production_and_oracle_are_distinct(self):
        for oracle in (False,True):
            with tempfile.TemporaryDirectory() as t:
                p=Path(t);self.fixture(p,oracle);r=evidence(p,2,True,oracle,1)
                self.assertEqual(r['production_retirement'],not oracle);self.assertEqual(r['zero_domain_ranks'],1)
    def test_bare_rejects_reference_output(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.assertEqual(evidence(p,2,False),{'enabled':False});self.fixture(p)
            with self.assertRaises(ValueError):evidence(p,2,False)
    def test_missing_transaction_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);f=p/'transform279-production-rank-0.txt';f.write_text('\n'.join(f.read_text().splitlines()[:-1]))
            with self.assertRaises(ValueError):evidence(p,2,True,False,1)
    def test_audit_tampering_rejected(self):
        for kind in ('reference','context','traffic','memory','work','epoch'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as t:
                p=Path(t);self.fixture(p)
                if kind=='work':
                    f=p/'transform279-work-rank-0.txt';f.write_text('1 1 1 2 3 4 6\n')
                elif kind=='epoch':
                    f=p/'transform279-cost-rank-0.txt';f.write_text(f.read_text()*2)
                else:
                    f=p/'transform279-production-rank-0.txt';rows=[x.split() for x in f.read_text().splitlines()]
                    index={'reference':8,'context':15,'traffic':9,'memory':13}[kind];rows[0][index]=str(int(rows[0][index])+1)
                    f.write_text('\n'.join(' '.join(x) for x in rows))
                with self.assertRaises(ValueError):evidence(p,2,True,False,1)
    def set_classes(self,root,rank,classes):
        f=root/f'transform279-production-rank-{rank}.txt'
        rows=[list(map(int,line.split())) for line in f.read_text().splitlines()]
        for row in rows:row[5:8]=classes
        f.write_text('\n'.join(' '.join(map(str,row)) for row in rows))
        (root/f'transform279-work-rank-{rank}.txt').write_text('1 1 1 2 '+' '.join(map(str,classes))+'\n')
    def test_outer_only_zero_domain_owner_is_valid(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);self.set_classes(p,1,[0,0,5]);r=evidence(p,2,True,False,1)
            self.assertEqual(r['zero_domain_kernel_records'],180)
            self.assertEqual(r['global_kernel_records_by_transaction'],[[3,4,10]]*36)
    def test_canonical_owner_with_no_records_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p);self.set_classes(p,1,[0,0,0])
            with self.assertRaisesRegex(ValueError,'no transform records'):evidence(p,2,True,False,1)
    def test_global_missing_class_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);self.fixture(p)
            for rank in range(2):self.set_classes(p,rank,[0,0,5])
            with self.assertRaisesRegex(ValueError,'Missing global'):evidence(p,2,True,False,1)
    def test_global_census_must_match_across_counts(self):
        a=dict(enabled=True,global_kernel_records_by_transaction=[[3,4,5]]*36)
        b=dict(enabled=True,global_kernel_records_by_transaction=[[3,4,6]]*36)
        check_global_census([a,dict(enabled=False),a])
        with self.assertRaisesRegex(ValueError,'census differs'):check_global_census([a,b])
    def test_selector_isolation(self):
        for kind in ('main','stage266','previous','stage279','candidate','checked'):
            e=environment(kind)
            self.assertEqual(e['WAVETRISK_RESIDENT_TRANSFORM'],str(int(kind in ('candidate','checked','stage279'))))
            self.assertEqual(e['WAVETRISK_RESIDENT_TRANSFORM_CHECK'],'0')
    def test_new_profile_regions_are_zero_for_predecessor(self):
        self.assertEqual(region([{'attribution':{'regions':{}}}],90),0)

if __name__=='__main__':unittest.main()
