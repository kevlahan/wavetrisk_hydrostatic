import struct
import tempfile
import unittest
from pathlib import Path
import numpy as np
from compare_physics_probe import load,compare


class PhysicsProbeTests(unittest.TestCase):
    def test_observed_rounding_midpoint_not_tolerance_waiver(self):
        legacy=0.06222941167660471
        block=0.062229411676523724
        packed=np.array([legacy,block],dtype='<f4')
        self.assertEqual(packed.view('u4').tolist(),[1031726149,1031726148])
        midpoint=float(packed.astype('f8').mean())
        self.assertLess(block,midpoint)
        self.assertGreater(legacy,midpoint)
        self.assertLess(legacy-block,1e-13)

    def record(self,value=1.):
        # Two atmospheric layers and one soil layer.
        h=struct.pack('<11i',181402,260,77,16,2,0,0,0,2,2,1)
        raw=np.arange(19,dtype='<f8').tobytes()
        packed=np.ones(23,dtype='<f4');packed[13]=value  # U(1)
        return h+raw+packed.tobytes()

    def test_fields_bits_and_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a';b=Path(tmp)/'b'
            a.write_bytes(self.record());b.write_bytes(self.record(np.nextafter(np.float32(1),np.float32(2))))
            self.assertEqual(compare(str(a),str(a))['fields']['U']['different'],0)
            c=compare(str(a),str(b))
            self.assertEqual(c['fields']['U']['different'],1)
            first=c['fields']['U']['first']
            self.assertEqual(first['index'],0)
            self.assertEqual(first['right_bits']-first['left_bits'],1)

    def test_reject_duplicate_truncated_nonfinite(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=Path(tmp)/'a'
            a.write_bytes(self.record()*2)
            with self.assertRaises(ValueError):load(str(a))
            a.write_bytes(self.record()[:-1])
            with self.assertRaises(ValueError):load(str(a))
            a.write_bytes(self.record(np.nan))
            with self.assertRaises(ValueError):compare(str(a),str(a))


if __name__=='__main__':unittest.main()
