import tempfile,unittest
from pathlib import Path
from contract import residency,footprint,READS
from test_contract import TEXT
from ownership import read_plan

class Inputs(unittest.TestCase):
    def test_masks_match_kernel_and_fortran_probe(self):
        import re,tarfile
        here=Path(__file__).resolve().parent
        path=here.parents[2]/'src/parallel_block_velocity.f90'
        if (here/'sources.tar.gz').is_file():
            with tarfile.open(here/'sources.tar.gz') as archive:
                text=archive.extractfile('candidate/src/parallel_block_velocity.f90').read().decode()
        else:text=path.read_text()
        body=text.split('  pure function velocity_qperp(',1)[1].split('  end function velocity_qperp',1)[0]
        body=body[body.index('value(1) ='):]
        for field in ('flux','pv'):
            reads={tuple(map(int,m)) for m in re.findall(r'\b'+field+r'\((\d+),(\d+)\)',body)}
            self.assertEqual(set(READS),reads)
        overlay=(here/'instrument.py').read_text()
        values=[]
        for name in ('component','position'):
            match=re.search(r'audit275_'+name+r'\(19\)=\[([^]]+)\]',overlay)
            values.append(tuple(map(int,match.group(1).split(','))))
        self.assertEqual(tuple(zip(*values)),READS)

    def fixture(self,root,change=None):
        for r in (0,1):
            rows=[f'RANK {r} 2 1 2']
            for k in (1,2):
                values=[10,20,30,40,50,60]
                if r==1 and change=='stale':values[2]=31
                if change!='missing' or r==0:
                    rows.append(('D' if r==0 else 'B')+f' 1 {k} 0 1 100 200 300 '+' '.join(map(str,values)))
                # Keep a second resident node so "missing" is a real incomplete match.
                if r==1:rows.append(f'B 1 {k} 0 1 101 200 300 10 20 30 40 50 60')
                rows.append(f'CHECK 1 {k} {k if r==0 else 0} {k if r==0 else 0}')
            (root/f'inputs275-rank-{r}.txt').write_text('\n'.join(rows)+'\n')
    def test_match_across_ranks_and_zero_domain_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root);r=residency(root,2)
            self.assertTrue(r['sampled_owner_interior_reuse_matches'])
            self.assertEqual(r['source_shadow_calls'],2)
            self.assertEqual(r['epochs'][0]['cross_rank_matches'],1)
    def test_stale_values_are_reported_not_credited(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root,'stale');r=residency(root,2)
            self.assertFalse(r['sampled_owner_interior_reuse_matches'])
            self.assertEqual(r['epochs'][0]['different_values']['rho'],1)
    def test_missing_values_are_reported_not_credited(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root,'missing');r=residency(root,2)
            self.assertFalse(r['sampled_owner_interior_reuse_matches'])
            self.assertEqual(r['epochs'][0]['missing_nodes'],1)
    def test_duplicate_coordinate_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root);p=root/'inputs275-rank-1.txt'
            p.write_text(p.read_text()+p.read_text().splitlines()[1]+'\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):residency(root,2)
    def test_truncated_rank_sequence_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.fixture(root);p=root/'inputs275-rank-1.txt'
            p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n')
            with self.assertRaises(ValueError):residency(root,2)
    def test_component_budget_deduplicates_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'plan';p.write_text(TEXT);s=read_plan(p)['snapshots'][0]
            f=footprint(s,0)
            self.assertEqual(f['flux']['logical_values'],6)
            self.assertEqual(f['flux']['Domain_remote_bytes'],24)
            self.assertEqual(f['physics']['logical_values'],6)
            self.assertEqual(f['rho']['unmapped_values'],1)
    def test_inactive_stencils_have_no_input_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'plan';p.write_text(TEXT);s=read_plan(p)['snapshots'][0]
            for a in s['actions']:
                if a['kind']=='DIRECT':a['active']=0
            for a in s['gradients']:a['active']=0
            self.assertEqual(footprint(s,0),{})

if __name__=='__main__':unittest.main()
