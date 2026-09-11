"""Keep both post-restart remaps mandatory in the arithmetic diagnostic."""
from pathlib import Path
import tempfile
import unittest

from run_arithmetic_diagnostic import coverage


class ArithmeticCoverage(unittest.TestCase):
    def check_log(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'run.log').write_text(text)
            return coverage(root)

    def test_both_restarts_have_following_remaps(self):
        text = ('Restarting from checkpoint 4\nRestarting from checkpoint 5\n'
                'Remapping vertical coordinates\nRestarting from checkpoint 6\n'
                'Remapping vertical coordinates\n')
        text += '00000322 0.33 d dt = 207.1 s\n' * 11
        self.assertEqual(self.check_log(text), {
            'reloads': [4, 5, 6], 'remap_after_each_new_restart': True, 'steps': 11})

    def test_remap_cannot_cover_the_wrong_restart(self):
        for text in [
            'Restarting from checkpoint 4\nRemapping vertical coordinates\n'
            'Restarting from checkpoint 5\nRestarting from checkpoint 6\n'
            'Remapping vertical coordinates\n',
            'Restarting from checkpoint 4\nRestarting from checkpoint 5\n'
            'Remapping vertical coordinates\nRestarting from checkpoint 6\n',
        ]:
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, 'Missing remap'):
                self.check_log(text)

    def test_extra_restart_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Expected checkpoint reloads'):
            self.check_log('Restarting from checkpoint 4\nRestarting from checkpoint 5\n'
                           'Restarting from checkpoint 5\nRestarting from checkpoint 6\n')


if __name__ == '__main__':
    unittest.main()
