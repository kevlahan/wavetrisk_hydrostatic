import json
from pathlib import Path
import tempfile
import unittest
from capability_223 import events, schedule
from placement_223 import plan
from collect_cluster_logs_223 import collect, verify_archive


class CapabilityTests(unittest.TestCase):
    def test_physical_cores_and_two_socket_mapping(self):
        rows=[]
        for s in range(2):
            for c in range(128):
                for thread in range(2):
                    rows.append(dict(cpu=s*1000+c*2+thread,core=c,socket=s,node=s*4+c//32,llc=[s*128+c//16]))
        cases=plan(dict(host='nk-node',rows=rows))
        for n,case in cases.items():
            mapping=case['mapping']
            self.assertEqual(len({r['cpu'] for r in mapping}),n)
            self.assertTrue(all(r['cpu']%2==0 for r in mapping))
            self.assertEqual([sum(r['socket']==s for r in mapping) for s in range(2)],[(n+1)//2,n//2])
        with self.assertRaises(ValueError):plan(dict(host='wrong',rows=rows[:256]))

    def test_gate_order(self):
        cases=list(schedule())
        self.assertEqual(len(cases),11)
        self.assertEqual([r[0] for r in cases],['J4']*6+['J5']*5)
        self.assertEqual(cases[0],('J4',40,'baseline',False,'native'))
        self.assertEqual(cases[5],('J4',41,'checked',False,'empty-root'))
        self.assertEqual([r[1] for r in cases[6:]],[160,160,167,167,167])

    def write_events(self, root):
        phases=('restart','install','dynamics','restart','install','step')
        for rank in range(41):
            lines=[]
            for seq,phase in enumerate(phases,1):
                lines.append(f'{phase} {seq} {rank} {int(rank!=0)} {-1 if phase=="restart" else 2} 1 3 100 1.00000000000000000E+00\n')
            (root/f'capability223-rank-{rank}.txt').write_text(''.join(lines))

    def test_empty_root_records_and_metadata_divergence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.write_events(root)
            e=events(root,41,40,1,'empty-root');self.assertEqual(e['restart_source_ranks'],[1])
            p=root/'capability223-rank-0.txt';p.write_text(p.read_text().replace(' 100 ',' 101 ',1))
            with self.assertRaisesRegex(ValueError,'metadata'):events(root,41,40,1,'empty-root')

    def test_no_idle_rank_and_missing_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.write_events(root)
            p=root/'capability223-rank-40.txt';p.write_text(p.read_text().replace('install 2 40 1 2','install 2 40 1 0'))
            with self.assertRaisesRegex(ValueError,'useful blocks'):events(root,41,40,1,'empty-root')
            p.unlink()
            with self.assertRaisesRegex(ValueError,'Incomplete'):events(root,41,40,1,'empty-root')

    def test_failed_run_archive_keeps_rank_evidence_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'results';out.mkdir();tree=root/'tree';tree.mkdir()
            (out/'results.json').write_text(json.dumps(dict(stage=223,status='failed',slurm_environment={'SLURM_JOB_ID':'123'})))
            (out/'capability223-rank-0.txt').write_text('partial evidence\n')
            archive=collect(223,'123',root,out,tree)
            record=verify_archive(archive)
            self.assertIn('results/capability223-rank-0.txt',record['files'])
            with self.assertRaises(FileExistsError):collect(223,'123',root,out,tree)


if __name__=='__main__':unittest.main()
