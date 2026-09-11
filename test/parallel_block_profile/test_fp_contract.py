"""Exercise the actual Makefile's arithmetic flags without compiling a solver."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class ArithmeticBuildPolicy(unittest.TestCase):
    @unittest.skipUnless(all(shutil.which(p) for p in ('make', 'gfortran', 'nf-config', 'nc-config')),
                         'GNU Fortran and NetCDF build tools required')
    def test_default_and_explicit_arithmetic_flags(self):
        makefile = Path(__file__).resolve().parents[2] / 'Makefile'
        env = {k: v for k, v in os.environ.items() if k not in ('FP_CONTRACT', 'FFLAGS', 'MAKEFLAGS')}
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / 'flags.mk'
            target.write_text("flags:\n\t@printf '%s\\n' '$(FFLAGS)'\n")
            for debug in ('false', 'check'):
                for choice, expected in ((None, 'off'), ('off', 'off'), ('default', None), ('fast', 'fast')):
                    with self.subTest(debug=debug, choice=choice):
                        command = ['make', '-s', '-f', str(makefile), '-f', str(target),
                                   'flags', 'PARAM=param_J5', 'DEBUG=' + debug]
                        if choice is not None:
                            command.append('FP_CONTRACT=' + choice)
                        result = subprocess.run(command, cwd=tmp, env=env, capture_output=True, text=True)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        flags = [s for s in result.stdout.split() if s.startswith('-ffp-contract=')]
                        self.assertEqual(flags, [] if expected is None else ['-ffp-contract=' + expected])
            result = subprocess.run(command + ['FP_CONTRACT=invalid'], cwd=tmp, env=env,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Unknown FP_CONTRACT', result.stderr)


if __name__ == '__main__':
    unittest.main()
