"""Compare the actual production polar callback with the legacy remap kernel."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


class PolarRemap(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'GNU Fortran required')
    def test_scalar_contract(self):
        helper = Path(__file__).resolve().parent
        source = (helper.parents[1]/'src/remap.f90').read_text()
        routines = []
        for name in ('remap_compressible_pole', 'remap_compressible', 'find_coordinates'):
            routines.append(re.search(r'^  subroutine '+name+r' \(.*?^  end subroutine '+name+r'$',
                                      source, re.M | re.S).group())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'polar_remap.inc').write_text('\n'.join(routines)+'\n')
            for optimization in ('-O0', '-O2'):
                subprocess.run(['gfortran', optimization, '-fcheck=all', '-ffree-line-length-132',
                                '-finit-real=snan', '-ffpe-trap=invalid,zero,overflow',
                                '-J'+tmp, '-I'+tmp, str(helper/'polar_remap.f90'),
                                '-o', str(root/'check')], cwd=root, check=True, capture_output=True)
                result = subprocess.run([str(root/'check')], check=True, capture_output=True, text=True)
                self.assertIn('polar remap contract passed', result.stdout)


if __name__ == '__main__':
    unittest.main()
