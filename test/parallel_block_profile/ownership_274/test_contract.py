import copy,tempfile,unittest
from pathlib import Path
from contract import read_plan,compile_snapshot,evidence

TEXT='''RANK 0 4 1 2
SNAP 0 1 7 4 1 2
NODE 1 0 0 2
NODE 2 1 1 2
NODE 3 1 1 1
NODE 4 -1 -1 -1
LEVEL 2 2
DIRECT 1 1 1 1 1 1 1 1 1 1 1
DIRECT 2 1 2 2 2 2 2 2 2 2 2
LEVEL 1 2
DIRECT 1 0 3 0 0 0 0 0 0 0 0
RESTRICT 2 3 1 2 0 0 1 0 0
GRADIENTS 2
GRADIENT 1 1 3 3 3 4
GRADIENT 2 0 4 0 0 0
END
'''

class Contract(unittest.TestCase):
    def parse(self,text=TEXT):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'plan.txt';p.write_text(text);return read_plan(p)
    def test_exact_dedup_budget_and_zero_domain_ranks(self):
        plan=self.parse();result=compile_snapshot(plan['snapshots'][0],0,4)
        self.assertEqual(result['source_input_budget']['logical_packed_values'],18)
        self.assertEqual(result['source_input_budget']['retained_Domain_producer_remote_bytes'],72)
        self.assertEqual(result['source_input_budget']['hypothetical_upstream_block_producer_remote_bytes'],0)
        self.assertEqual(result['restriction_remote_versioned_values'],1)
        self.assertEqual(result['same_level_cross_owner_dependencies'],0)
        self.assertEqual(result['unresolved_source_reads'],0)
        self.assertTrue(result['no_detected_same_level_or_uninitialized_source_blocker'])
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['unmapped_target_adapter_work'],{'gradient':1})
        self.assertEqual(result['gradient_input_budget']['unresolved_adapter_input_values'],4)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'ownership274-rank-0.txt').write_text(TEXT)
            for rank in range(1,4):(root/f'ownership274-rank-{rank}.txt').write_text(f'RANK {rank} 4 1 2\n')
            result=evidence(root,4,True,domains=1)
            self.assertEqual(result['by_sample']['1']['work_by_rank'][1]['direct'],2)
            self.assertFalse(result['performance_qualified'])
            with self.assertRaises(ValueError):evidence(root,4,False,domains=1)
            (root/'ownership274-rank-3.txt').unlink()
            with self.assertRaises(ValueError):evidence(root,4,True,domains=1)
    def test_same_level_exchange_is_a_design_blocker(self):
        snapshot=self.parse()['snapshots'][0]
        snapshot['actions'][-1]['level']=2
        result=compile_snapshot(snapshot,0,4)
        self.assertEqual(result['same_level_cross_owner_dependencies'],1)
        self.assertFalse(result['no_detected_same_level_or_uninitialized_source_blocker'])
    def test_overwrite_version_and_read_before_write(self):
        snapshot=self.parse()['snapshots'][0]
        # Same target is first directly zeroed, then partially restricted.
        # Its untouched components are defined by that direct write.
        snapshot['actions'].pop(2)
        result=compile_snapshot(snapshot,0,4)
        self.assertEqual(result['unresolved_source_reads'],2)
        self.assertFalse(result['no_detected_same_level_or_uninitialized_source_blocker'])
    def test_inactive_does_not_read_invalid_stencil(self):
        snapshot=self.parse()['snapshots'][0]
        snapshot['actions'][0]['active']=0
        result=compile_snapshot(snapshot,0,4)
        self.assertEqual(result['source_input_budget']['logical_packed_values'],9)
    def test_reject_corrupted_plan_records(self):
        variants=[TEXT.replace('NODE 2 1','NODE 1 1'),TEXT.replace('NODE 2 1','NODE 2 4'),
            TEXT.replace('NODE 4 -1 -1 -1','NODE 4 0 -1 -1'),TEXT.replace('LEVEL 2 2','LEVEL 3 2'),
            TEXT.replace('DIRECT 1 1 1','DIRECT 1 1 9'),TEXT.replace('RESTRICT 2','RESTRICT 3'),
            TEXT.replace('GRADIENT 1 1 3','GRADIENT 1 1 9'),TEXT.replace('SNAP 0 1','SNAP 0 2'),
            TEXT.removesuffix('END\n'),TEXT.replace(' 1 0 0\nGRADIENTS',' 2 0 0\nGRADIENTS')]
        for value in variants:
            with self.subTest(value=value[-80:]),self.assertRaises((ValueError,StopIteration)):self.parse(value)
    def test_rank_local_samples_survive_domain_migration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'ownership274-rank-0.txt').write_text(TEXT)
            (root/'ownership274-rank-1.txt').write_text(TEXT.replace('RANK 0','RANK 1'))
            for r in (2,3):(root/f'ownership274-rank-{r}.txt').write_text(f'RANK {r} 4 1 2\n')
            result=evidence(root,4,True,domains=1)
            self.assertEqual(result['domains_observed_on_multiple_ranks'],{'0':[0,1]})
            self.assertEqual(result['sampled_plans'],2)
            self.assertEqual(result['by_sample']['1']['domains'],1)
            self.assertEqual(result['by_sample']['1']['observations'],2)
            self.assertEqual([p['observing_rank'] for p in result['plans']],[0,1])
            self.assertNotEqual(result['plans'][0]['gradient_input_budget'],result['plans'][1]['gradient_input_budget'])
            # Duplicate local ordinals remain corruption, even at the same generation.
            (root/'ownership274-rank-1.txt').write_text(TEXT.replace('RANK 0','RANK 1')+'SNAP'+TEXT.split('SNAP',1)[1])
            with self.assertRaises(ValueError):evidence(root,4,True,domains=1)

if __name__=='__main__':unittest.main()
