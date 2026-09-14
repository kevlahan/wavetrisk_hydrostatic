import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from cluster_189 import finish_cross_mode


class Cluster189(unittest.TestCase):
    def test_cross_mode_gates_are_required_and_recorded(self):
        for failed in (None, 0, 1):
            with self.subTest(failed=failed), tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp)
                (out / 'results.json').write_text(json.dumps({
                    'status': 'passed', 'parameters': {'run_id': 'test'}, 'gates': {},
                    'limitations': ['old limitation']}))
                gates = [dict(passed=i != failed, failures=[] if i != failed else ['different']) for i in range(2)]
                with patch('cluster_189.exact_gate', side_effect=gates) as compare:
                    if failed is None:
                        finish_cross_mode(out)
                    else:
                        with self.assertRaisesRegex(ValueError, 'checkpoint mismatch'):
                            finish_cross_mode(out)
                    result = json.loads((out / 'results.json').read_text())
                    self.assertEqual(result['status'], 'passed' if failed is None else 'failed')
                    for call in compare.call_args_list:
                        self.assertEqual(call.args[0], out / 'checked-off')
                        self.assertEqual(call.args[2:], ('test', [4]))
                    self.assertEqual(len(result['gates']), 2 if failed is None else failed + 1)

    def test_does_not_accept_failed_campaign(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / 'results.json').write_text('{"status":"failed"}')
            with patch('cluster_189.exact_gate') as compare:
                with self.assertRaisesRegex(ValueError, 'did not pass'):
                    finish_cross_mode(out)
                compare.assert_not_called()
