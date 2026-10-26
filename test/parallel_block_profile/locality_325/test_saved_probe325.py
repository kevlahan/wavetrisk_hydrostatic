import tempfile
import unittest
from pathlib import Path
import probe_saved_matrix325 as p
import verify322 as v

class SavedProbeCensusTests(unittest.TestCase):
    def test_explicit_capture_and_report_paths_resolve(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);base,output=p.resolve_probe_paths(root/'captures',root/'review/report.json')
            self.assertEqual(base,(root/'captures').resolve())
            self.assertEqual(output,(root/'review/report.json').resolve())
    def test_report_cannot_overwrite_capture_tree(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for output in (root,root/'report.json',root/'qual-C0/execution316-rank-1.bin.gz'):
                with self.assertRaisesRegex(v.AuditError,'output_in_capture_root'):p.resolve_probe_paths(root,output)
    def fixture(self,root):
        for arm,ranks in p.EXPECTED.items():
            folder=root/arm;folder.mkdir()
            for rank in ranks:(folder/f'execution316-rank-{rank}.bin.gz').write_bytes(b'fixture')
    def test_absent_capture_root_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(v.AuditError,'exact_capture_inventory'):p.capture_inventory(Path(d)/'absent')
    def test_empty_capture_root_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(v.AuditError,'exact_capture_inventory'):p.capture_inventory(d)
    def test_exact_nine_capture_names_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            self.fixture(Path(d));self.assertEqual([(a,r) for a,r,path in p.capture_inventory(d)],[(a,r) for a,rs in p.EXPECTED.items() for r in rs])
    def test_missing_capture_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            self.fixture(Path(d));(Path(d)/'qual-C0/execution316-rank-2.bin.gz').unlink()
            with self.assertRaisesRegex(v.AuditError,'exact_capture_inventory'):p.capture_inventory(d)
    def test_extra_capture_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            self.fixture(Path(d));(Path(d)/'qual-C0/execution316-rank-3.bin.gz').write_bytes(b'fixture')
            with self.assertRaisesRegex(v.AuditError,'exact_capture_inventory'):p.capture_inventory(d)
    def test_empty_common_feature_set_rejected(self):
        for a,b in [({},{}),({'a':'x'},{'b':'x'})]:
            with self.assertRaisesRegex(v.AuditError,'empty_common_features'):p.common_feature_keys(a,b)

if __name__=='__main__':unittest.main()
