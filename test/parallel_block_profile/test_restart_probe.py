"""Synthetic layout and sensitivity tests for the read-only restart comparator."""
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compare_restart_probe import compare, load


class RestartProbeTests(unittest.TestCase):
    def write_probe(self, path, values, gid=7):
        # Two layers: k=0 (soil), k=1 (atmosphere), one active patch.
        raw = struct.pack('<8id', 179401, 1, 3, 0, 1, 4, 259, 1, 100.)
        raw += np.zeros(6, dtype='<f8').tobytes()
        raw += struct.pack('<2iqi', gid, 1, 1, 5)
        raw += np.full(64, 2, dtype='<i4').tobytes()
        raw += np.asarray(values, dtype='<f8').tobytes()
        path.write_bytes(raw)

    def test_identity_and_mass_field_location(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/'a.bin', Path(tmp)/'b.bin'
            values = np.zeros((2, 160))
            self.write_probe(a, values)
            same = compare(str(a), str(a))
            self.assertEqual(same['common_patches'], 1)
            self.assertTrue(all(s['different'] == 0 for s in same['fields'].values()))
            values[1, 99] = 3.5  # Atmospheric solution mass, cell index 3.
            self.write_probe(b, values)
            result = compare(str(a), str(b))
            mass = result['fields']['sol/mass/atmosphere/active']
            self.assertEqual(mass['different'], 1)
            self.assertEqual(mass['max_abs'], 3.5)
            self.assertEqual(mass['worst']['index'], 3)
            self.assertEqual(mass['worst']['k'], 1)
            self.assertEqual(result['fields']['wav/mass/atmosphere/all']['different'], 0)

    def test_missing_patch_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/'a.bin', Path(tmp)/'b.bin'
            values = np.zeros((2, 160))
            self.write_probe(a, values)
            self.write_probe(b, values, gid=8)
            result = compare(str(a), str(b))
            self.assertEqual(result['left_only_patches'], 1)
            self.assertEqual(result['right_only_patches'], 1)
            values[1, 128] = np.nan
            self.write_probe(b, values)
            result = compare(str(a), str(b))
            self.assertEqual(result['fields']['sol/temperature/atmosphere/active']['nonfinite'], 1)

    def test_no_files_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                load(str(Path(tmp)/'missing-*.bin'))


if __name__ == '__main__':
    unittest.main()
