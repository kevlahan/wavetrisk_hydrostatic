from pathlib import Path
import unittest
import io
import tarfile
import tempfile
import json
from unittest.mock import patch
import numpy as np
from local_j4 import refinement, vtk_refinement, compare_windows


class J4RefinementTests(unittest.TestCase):
    def fixture(self):
        threshold=np.tile([1.,2.,3.],41)
        waves=np.zeros((6,41,80));flags=np.zeros((6,4),dtype=int)
        flags[0,0]=1;flags[4,0]=1
        waves[5,11,0]=1.5;waves[5,11,1]=0.8
        waves[5,11,48]=2.5;waves[5,11,49]=1.5
        waves[5,11,64]=2.5;waves[5,11,65]=3.5
        return ((),threshold,np.zeros(2),np.zeros(80),waves,flags)

    def test_finest_threshold_counts_use_velocity_mass_temperature_order(self):
        with patch('local_j4.read',return_value=[self.fixture()]) as reader,patch('local_j4.digest',return_value='hash'):
            r=refinement(Path('test'))
        reader.assert_called_once_with(Path('test'),-10,30,expected_domains=40)
        self.assertEqual(r['patches_by_level'],{4:4,5:1,6:1})
        self.assertEqual(r['above_threshold_node_wavelet_locations'][6],2)
        self.assertEqual(r['above_threshold_edge_wavelet_locations'][6],1)
        self.assertEqual(r['finest_above_threshold_locations'],3)
        self.assertEqual(r['finest_fraction_of_above_threshold_locations'],1.)

    def test_missing_tree_and_nonfinite_rejected(self):
        record=self.fixture();record[-1][5,0]=1
        with patch('local_j4.read',return_value=[record]),patch('local_j4.digest',return_value='hash'):
            with self.assertRaisesRegex(ValueError,'Missing'):refinement(Path('test'))
        record=self.fixture();record[-2][0,11,0]=np.nan
        with patch('local_j4.read',return_value=[record]):
            with self.assertRaisesRegex(ValueError,'Nonfinite'):refinement(Path('test'))

    def test_exported_active_mesh_levels_and_vertices(self):
        header=b'# vtk DataFile Version 2.0\nWAVETRISK adaptive data\nBINARY\nDATASET POLYDATA\nPOINTS 4 float\n'
        raw=header+np.zeros((4,3),dtype='>f4').tobytes()+b'POLYGONS 2 8\n'
        raw+=np.array([[3,0,1,2],[3,1,2,3]],dtype='>i4').tobytes()
        raw+=b'CELL_DATA 2\nSCALARS Level int\nLOOKUP_TABLE default\n'+np.array([4,6],dtype='>i4').tobytes()
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mesh.tgz'
            with tarfile.open(path,'w:gz') as archive:
                info=tarfile.TarInfo('test_tri_001_0003.vtk');info.size=len(raw)
                archive.addfile(info,io.BytesIO(raw))
            r=vtk_refinement(path)
            self.assertEqual(r['active_triangles_by_level'],{4:1,6:1})
            self.assertEqual(r['vertices_on_level6_triangles'],3)
            self.assertEqual(r['level6_triangle_fraction'],0.5)

    def test_single_restart_window_keeps_exact_gate_and_remap_requirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            left=Path(tmp)/'a';right=Path(tmp)/'b'
            identity={'parameters':{'run_id':'test'},'input_sha256':'input','seed_sha256':'seed','grids':{},'ranks':4}
            for root in (left,right):
                root.mkdir();(root/'results.json').write_text(json.dumps(identity))
                (root/'test_checkpoint_0004.bin.zst').touch()
                (root/'run.log').write_text('Restarting from checkpoint 4\nRemapping vertical coordinates\nTotal cpu time = 1.0\n')
            report={'header_differences':0,'topology_differences':0,'threshold_nonfinite':0,'threshold_max_abs':0.,
                    'fields':{'coarse/mass/atmosphere':{'max_abs':0.,'nonfinite':0}}}
            with patch('local_j4.compare',return_value=report) as compare:
                self.assertTrue(compare_windows(left,right,'seed',expected_new=(4,))['passed'])
                self.assertEqual(compare.call_args.kwargs,{'expected_domains':40})
                report['fields']['coarse/mass/atmosphere']['max_abs']=1e-15
                self.assertFalse(compare_windows(left,right,'seed',expected_new=(4,))['passed'])
                report['fields']['coarse/mass/atmosphere']['max_abs']=0.
                (right/'run.log').write_text('Restarting from checkpoint 4\nTotal cpu time = 1.0\n')
                self.assertFalse(compare_windows(left,right,'seed',expected_new=(4,))['passed'])
            identity['input_sha256']='different';(right/'results.json').write_text(json.dumps(identity))
            with self.assertRaisesRegex(ValueError,'Unmatched'):compare_windows(left,right,'seed',expected_new=(4,))


if __name__=='__main__':unittest.main()
