"""Check independent post-MPI verification detects the observed failure classes."""
from pathlib import Path
import struct
import tempfile
import unittest
from verify_checkpoint_io_probe import verify


class ProbeVerifierTests(unittest.TestCase):
    def test_valid_zero_header_and_wrong_rank_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'probe.bin'
            header = struct.pack('<3q', 0x5741564554524953, 1, 160)
            directory = (struct.pack('<160i', *([1] * 160)) +
                         struct.pack('<160q', *(3224 + 16 * i for i in range(160))) +
                         struct.pack('<160q', *([16] * 160)))
            original = bytearray(header + directory + bytes([1]) * 8192 + bytes([2]) * 8192)
            p.write_bytes(original)
            r = verify(p, 2)
            self.assertFalse(r['bad_header'] or r['bad_directory'] or r['bad_payload_blocks'] or r['bad_size'])
            damaged = original.copy()
            damaged[:24] = bytes(24)
            damaged[3224 + 8192 + 17] = 0
            p.write_bytes(damaged)
            r = verify(p, 2)
            self.assertTrue(r['bad_header'])
            self.assertFalse(r['bad_directory'])
            self.assertEqual(r['bad_payload_blocks'], [{'rank_block': 1, 'file_offset': 3224 + 8192 + 17,
                             'mismatched_or_missing_bytes': 1, 'expected_byte': 2, 'observed_byte': 0}])
            p.write_bytes(original[:-19])
            r = verify(p, 2)
            self.assertTrue(r['bad_size'])
            self.assertEqual(r['bad_payload_blocks'][0]['mismatched_or_missing_bytes'], 19)
            self.assertIsNone(r['bad_payload_blocks'][0]['observed_byte'])


if __name__ == '__main__':
    unittest.main()
