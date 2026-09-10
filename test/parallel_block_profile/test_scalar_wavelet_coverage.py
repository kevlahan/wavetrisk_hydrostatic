"""Compile the actual producer with minimal stubs to test traversal coverage."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


class ScalarWaveletCoverage(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'GNU Fortran required')
    def test_production_traversal(self):
        helper = Path(__file__).resolve().parent
        source = (helper.parents[1]/'src/parallel_block_mpi.f90').read_text()
        producer = re.search(r'^  subroutine produce_block_scalar_wavelets .*?^  end subroutine produce_block_scalar_wavelets',
                             source, re.M | re.S).group()
        context = re.search(r'^  type :: Block_Scalar_Wavelet_Context.*?^  end type Block_Scalar_Wavelet_Context',
                            source, re.M | re.S).group()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'producer.inc').write_text(producer+'\n')
            (root/'context.inc').write_text(context+'\n')
            for optimization in ('-O0', '-O2'):
                subprocess.run(['gfortran', optimization, '-fcheck=all', '-ffree-line-length-132',
                                '-J'+tmp, '-I'+tmp, str(helper/'scalar_wavelet_coverage.f90'),
                                '-o', str(root/'check')], cwd=root, check=True, capture_output=True)
                result = subprocess.run([str(root/'check')], check=True, capture_output=True, text=True)
                self.assertIn('coverage passed', result.stdout)
            # A producer omitting the external-parent root must fail the
            # regression, even though its visited sites are internally valid.
            broken = producer.replace('do p = 1,size(block%patch)', 'do p = 2,size(block%patch)', 1)
            self.assertNotEqual(broken, producer)
            (root/'producer.inc').write_text(broken+'\n')
            subprocess.run(['gfortran', '-O0', '-fcheck=all', '-J'+tmp, '-I'+tmp,
                            str(helper/'scalar_wavelet_coverage.f90'), '-o', str(root/'broken')],
                           cwd=root, check=True, capture_output=True)
            result = subprocess.run([str(root/'broken')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('missing root coverage', result.stderr)


if __name__ == '__main__':
    unittest.main()
