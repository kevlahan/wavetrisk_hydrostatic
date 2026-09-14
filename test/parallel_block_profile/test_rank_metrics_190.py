import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from rank_metrics_190 import rank_resources,snapshot_allocation


class RankMetrics190(unittest.TestCase):
    def test_resource_units_and_missing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'rank-0.json').write_text('{"rank":0}')
            with self.assertRaises(FileNotFoundError):rank_resources(root)
            (root/'rank-0.resources').write_text('Maximum resident set size (kbytes): 2048\n'
                'Major (requiring I/O) page faults: 3\nMinor (reclaiming a frame) page faults: 7\n')
            row=rank_resources(root)[0]
            self.assertEqual((row['max_rss_kib'],row['major_faults'],row['minor_faults']),(2048,3,7))

    def test_overlapping_rank_affinity_is_rejected(self):
        rows=[dict(rank=i,host='node',allowed_cpus=[0]) for i in range(2)]
        with tempfile.TemporaryDirectory() as tmp,patch('rank_metrics_190.collect',return_value=rows), \
             patch('rank_metrics_190.subprocess.run',return_value=subprocess.CompletedProcess([],0,stdout='CPU topology')):
            with self.assertRaisesRegex(ValueError,'Overlapping'):snapshot_allocation(Path(tmp),2)
            rows[1]['allowed_cpus']=[1]
            self.assertEqual(snapshot_allocation(Path(tmp),2)['ranks'],rows)

if __name__=='__main__':unittest.main()
