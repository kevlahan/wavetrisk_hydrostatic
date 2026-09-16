"""Exercise extracted production interpolation against the original scalar path."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def extract(source, kind, name):
    return re.search(r'^  (?:real\(dp\) )?' + kind + r' ' + name +
                     r'\s*\(.*?^  end ' + kind + ' ' + name + r'\b', source, re.M | re.S).group()


def prepare(out):
    source = (ROOT / 'src/parallel_block_mpi.f90').read_text()
    shared = (ROOT / 'src/shared.f90').read_text()
    names = ('interpolate_block_scalar_flux', 'interpolate_block_scalar_flux_columns',
             'block_restriction_edge_flux')
    for name in names:
        (out / (name + '.inc')).write_text(extract(source, 'function', name) + '\n')
    (out / 'indices.inc').write_text(extract(source, 'subroutine', 'get_block_restriction_indices') + '\n')
    constants = []
    for name in ('UZM UPZ UMZ UZP VMM VPM VMP VPP WMM WPM WMP WPP WMMM WPPP VMPP VPMM MP PP PM MM '
                 '').split():
        constants.append(re.search(r'^  integer, parameter :: ' + name + r'\s*=.*$', shared, re.M).group())
    for name in ('RT', 'EDGE'):
        constants.append(re.search(r'^  integer, parameter :: .*\b' + name + r'\s*=.*$', shared, re.M).group())
    block = (ROOT / 'src/parallel_block.f90').read_text()
    for name in ('STORE_PATCH', 'STORE_BDRY', 'STORE_GHOST'):
        constants.append(re.search(r'^  integer, parameter, public :: ' + name + r'\s*=.*$', block, re.M).group())
    constants.append(re.search(r'^  real\(dp\), parameter :: BLOCK_RESTRICTION_VALUE_LIMIT = &\n[^\n]*',
                               source, re.M).group())
    for name in ('hex_sides', 'hex_s_offs', 'nghb_pt'):
        lines = shared.splitlines()
        index = next(i for i, line in enumerate(lines) if re.match(
            r'  integer, parameter :: ' + name + r'\b', line))
        declaration = [lines[index]]
        while declaration[-1].rstrip().endswith('&'):
            index += 1
            declaration.append(lines[index])
        constants.append('\n'.join(declaration))
    for name in ('BLOCK_SCALAR_DIRECT_FLUX_START', 'BLOCK_SCALAR_RESTRICTION_WEIGHT_START'):
        constants.append(re.search(r'^  integer, parameter :: ' + name + r'\s*=.*$', source, re.M).group())
    (out / 'constants.inc').write_text('\n'.join(constants) + '\n')
    shutil.copy2(Path(__file__).with_name('restriction_columns.f90'), out)


class RestrictionColumns(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'GNU Fortran required')
    def test_interpolation_and_invalid_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            prepare(out)
            for opt in ('-O0', '-O2'):
                build = subprocess.run(['gfortran','-cpp', opt, '-ffp-contract=off', '-fcheck=all',
                                        '-ffpe-trap=invalid,zero,overflow', '-Werror=line-truncation',
                                        '-I' + tmp, str(ROOT / 'src/kind.f90'),
                                        str(ROOT / 'src/parallel_block_scalar_storage.f90'),
                                        'restriction_columns.f90', '-o', 'check'], cwd=out,
                                       text=True, capture_output=True)
                self.assertEqual(build.returncode, 0, build.stdout + build.stderr)
                for mode in ('valid', 'poison', 'nan', 'extent'):
                    run = subprocess.run(['./check', mode], cwd=out, text=True, capture_output=True)
                    if mode == 'valid':
                        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                        self.assertIn('PASS restriction columns', run.stdout)
                    else:
                        self.assertNotEqual(run.returncode, 0, mode)
                        self.assertIn({'poison': 'column is poisoned', 'nan': 'column is non-finite',
                                       'extent': 'column'}[mode], run.stdout + run.stderr)


if __name__ == '__main__':
    unittest.main()
