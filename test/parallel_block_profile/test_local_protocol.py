import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from analyze import analyze
from checkpoint_gate import accept, topology, compare_runs, POLICY
from local_campaign import require_certificate
from profile_protocol import child_rss, environment, grid_identity, checkpoint_name, parameters
from prepare_legacy_profile import boundary_patch
from summarize_local_campaign import aggregate


class ProtocolTests(unittest.TestCase):
    def test_nested_grid_links_hash_contents_and_reject_cycles(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);grid=root/'grid';asset=root/'assets'
            grid.mkdir();asset.mkdir();(asset/'a').write_text('first')
            (grid/'linked').symlink_to(asset,target_is_directory=True)
            first=grid_identity(grid)
            self.assertEqual(first['files'],1)
            (asset/'a').write_text('second')
            self.assertNotEqual(first,grid_identity(grid))
            (asset/'cycle').symlink_to(grid,target_is_directory=True)
            with self.assertRaisesRegex(ValueError,'cycle'):grid_identity(grid)

    def test_dangling_grid_link_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'broken').symlink_to(root/'absent')
            with self.assertRaises(FileNotFoundError):grid_identity(root)

    def test_environment_resets_inherited_switches(self):
        with patch.dict('os.environ',{'WAVETRISK_UNKNOWN':'1','WAVETRISK_VALIDATE_BLOCK_REMAP':'1'}):
            off=environment();on=environment(oracles=True)
        self.assertNotIn('WAVETRISK_UNKNOWN',off)
        for flag in ('DYNAMICS','ADAPTATION','REMAP'):
            self.assertEqual(off['WAVETRISK_VALIDATE_BLOCK_'+flag],'0')
            self.assertEqual(on['WAVETRISK_VALIDATE_BLOCK_'+flag],'1')
        self.assertEqual(off['OMP_NUM_THREADS'],'1')

    def test_rss_tracks_only_launcher_descendants(self):
        rows='10 1 20 mpirun\n11 10 40 daemon\n12 11 600 /tmp/climateJ5\n13 1 900 /tmp/climateJ5\n'
        self.assertEqual(child_rss(rows,10),{'12':600})

    def test_seed_is_not_relabelled(self):
        self.assertEqual(checkpoint_name('run_id test\nresume 4\n'),'test_checkpoint_0004.bin.zst')
        with self.assertRaises(ValueError):parameters('resume 3\nresume 4\n')
        with self.assertRaises(ValueError):checkpoint_name('run_id ../escape\nresume 4\n')

    def test_remap_and_checkpoint_are_separate_step_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'run.log'
            path.write_text('Remapping vertical coordinates ...\n'
                            '00000322 0.3 d dt = 1 s cpu = 2.0\n'
                            'Saving checkpoint 4\n00000322 0.4 d dt = 1 s cpu = 3.0\n'
                            '00000322 0.5 d dt = 1 s cpu = 1.0\nTotal cpu time = 6.0\n')
            rows=analyze(path)['step_records']
        self.assertEqual([(r['checkpoint'],r['remap']) for r in rows],[(False,True),(True,False),(False,False)])

    def test_boundary_patch_is_pinned_and_nonempty(self):
        repo=Path(__file__).resolve().parents[2]
        data=boundary_patch(repo)
        self.assertIn(b'parallel_block_profile_mod',data)

    def test_critical_rank_attribution_conserves_mean_step(self):
        a=[0.0]*65;b=[0.0]*65
        a[0]=1;a[1]=3;b[0]=5;b[1]=1
        summary=aggregate([{'wall':4,'cpu':2,'self':a},{'wall':6,'cpu':3,'self':b}],{})
        self.assertEqual(summary['mean_max_rank_wall'],5)
        self.assertEqual(sum(r['self_wall_on_step_critical_rank'] for r in summary['regions']),5)


class GateTests(unittest.TestCase):
    def report(self):
        return {'header_differences':0,'topology_differences':0,'threshold_nonfinite':0,
                'threshold_max_abs':0,'fields':{'pole/mass/atmosphere':{'nonfinite':0,'max_abs':1e-12}}}

    def test_finite_small_differences_are_not_claimed_exact(self):
        self.assertFalse(accept(self.report()))
        self.assertTrue(accept(self.report(),exact=True))

    def test_pole_error_and_nonfinite_threshold_fail(self):
        report=self.report();report['fields']['pole/mass/atmosphere']['max_abs']=0.1
        self.assertTrue(accept(report))
        report=self.report();report['threshold_nonfinite']=1
        self.assertTrue(accept(report))

    def test_soil_and_topology_must_agree(self):
        report=self.report();report['fields']={'coarse/temperature/soil':{'nonfinite':0,'max_abs':1e-12}}
        self.assertTrue(accept(report))
        report=self.report();report['topology_differences']=1
        self.assertTrue(accept(report))

    def test_topology_coverage_and_nonfinite(self):
        waves=np.ones((4,41,80));flags=np.zeros((4,4),dtype=int)
        record=((),np.zeros(3),np.zeros(2),np.zeros(80),waves,flags)
        with patch('checkpoint_gate.read',return_value=[record]):
            self.assertEqual(topology(Path('unused'))['maximum_level'],5)
        flags[0,0]=1
        with patch('checkpoint_gate.read',return_value=[record]):
            with self.assertRaisesRegex(ValueError,'Missing'):topology(Path('unused'))
        waves[0,0,0]=np.nan
        with patch('checkpoint_gate.read',return_value=[record]):
            with self.assertRaisesRegex(ValueError,'Nonfinite'):topology(Path('unused'))

    def test_certificate_is_bound_and_fail_closed(self):
        identity={'binaries':{'block':'abc'},'fixture':{'parameters':{'resume':'4'}},'ranks':4}
        cert={'passed':True,'policy':POLICY,'binding':copy.deepcopy(identity),'validation_time_end':'0.335'}
        require_certificate(cert,identity,'0.319')
        with self.assertRaises(ValueError):require_certificate({},identity,'0.319')
        with self.assertRaises(ValueError):require_certificate(cert,identity,'0.34')
        identity['ranks']=8
        with self.assertRaises(ValueError):require_certificate(cert,identity,'0.319')
        identity['ranks']=4;identity['binaries']['block']='changed'
        with self.assertRaises(ValueError):require_certificate(cert,identity,'0.319')

    def test_checkpoint_cycle_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,'at least two'):compare_runs(root,root,'seed.zst')
            for n in (5,6):(root/f'test_checkpoint_{n:04d}.bin.zst').touch()
            (root/'run.log').write_text('Restarting from checkpoint 5\n')
            with patch('checkpoint_gate.topology',return_value={'patches_by_level':{5:640}}), \
                 patch('checkpoint_gate.compare',return_value=self.report()):
                result=compare_runs(root,root,'seed.zst')
            self.assertFalse(result['passed'])
            self.assertTrue(any('no remap' in reason for reason in result['failures']))
            self.assertTrue(any('missing reload 6' in reason for reason in result['failures']))


if __name__=='__main__':unittest.main()
