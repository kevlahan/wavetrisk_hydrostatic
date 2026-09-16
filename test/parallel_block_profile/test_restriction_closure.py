"""Exact extracted closure comparison with real storage and synthetic routing."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from test_restriction_columns import ROOT, extract, prepare as prepare_interpolation


def prepare(out):
    prepare_interpolation(out)
    source = (ROOT / 'src/parallel_block_mpi.f90').read_text()
    shared = (ROOT / 'src/shared.f90').read_text()
    constants = (out / 'constants.inc').read_text()
    for name in ('BLOCK_SCALAR_NATIVE_DSCALAR_INDEX', 'BLOCK_SCALAR_RESTRICTED_FLUX_START',
                 'BLOCK_SCALAR_AREA_INDEX', 'BLOCK_SCALAR_ACTIVE_INDEX', 'BLOCK_SCALAR_OVERLAP_START',
                 'BLOCK_SCALAR_FLUX_COUNT'):
        constants += re.search(r'^  integer, parameter :: ' + name + r'\s*=.*$', source, re.M).group() + '\n'
    lines=shared.splitlines()
    i=next(i for i,l in enumerate(lines) if re.match(r'  integer, parameter :: bfly_no2\b',l))
    decl=[lines[i]]
    while decl[-1].rstrip().endswith('&'):
        i+=1; decl.append(lines[i])
    (out/'constants.inc').write_text(constants+'\n'.join(decl)+'\n')
    names=('complete_block_coarse_flux', 'partial_block_coarse_flux', 'block_coarse_flux',
           'block_area_value', 'block_overlap_value', 'block_overlap_split', 'block_scalar_flux_divergence',
           'complete_block_coarse_flux_columns')
    for name in names:
        (out/(name+'.inc')).write_text(extract(source,'function',name)+'\n')
    # Unchanged scalar production core; remove the observational oracle branch
    # whose Domain/MPI provenance is covered by solver oracle runs, not this fixture.
    ds=extract(source,'function','block_restriction_dscalar')
    ds=ds[:ds.index('    ! Stage 166B')]+ '  end function block_restriction_dscalar\n'
    (out/'block_restriction_dscalar.inc').write_text(ds)
    for name in ('validate_scalar_restriction_vector','validate_scalar_restriction_value'):
        (out/(name+'.inc')).write_text(extract(source,'subroutine',name)+'\n')
    original=(out/'restriction_columns.f90').read_text().split('program check')[0]
    original=original.replace('only: int64','only: int64,error_unit')
    original=original.replace('  type :: Block_Data','  type :: Block_Scalar_Restriction_Cursor\n     integer :: unused\n  end type\n  type :: Block_Data')
    original=original.replace('inverse_address_work(8)','inverse_address_work(12)')
    original=original.replace('x,y,slot) result(value)','x,y,slot,allow_uninitialized) result(value)')
    original=original.replace('    real(dp) :: value\n','    logical, optional, intent(in) :: allow_uninitialized\n    logical :: permit\n    real(dp) :: value\n')
    original=original.replace('    scalar_reads=scalar_reads+1','''    permit=.false.
    if(present(allow_uninitialized))permit=allow_uninitialized
    if(.not.ieee_is_finite(value))call fail('scalar value is non-finite')
    if(.not.permit.and.abs(value)>=BLOCK_RESTRICTION_VALUE_LIMIT)call fail('scalar value is poisoned')
    scalar_reads=scalar_reads+1''')
    original=original.replace('       call scalar_write_range', '''       values(7)=0.25_dp+real(node+tag,dp)/901.0_dp
       values(8)=1.0_dp
       values(21:26)=0.1_dp+real(node+tag,dp)/997.0_dp
       ! Deliberately mix valid transported tendencies and missing layers.
       if(modulo(layer,3)==0)values(BLOCK_SCALAR_NATIVE_DSCALAR_INDEX)=-0.125_dp*huge(0.0_dp)
       call scalar_write_range''')
    includes='\n'.join("  include '"+n+".inc'" for n in (*names,'block_restriction_dscalar',
        'validate_scalar_restriction_vector','validate_scalar_restriction_value'))
    original=original.replace('end module', '''  subroutine report_scalar_restriction_cursor(message,cursor)
    character(*),intent(in) :: message
    type(Block_Scalar_Restriction_Cursor),intent(in) :: cursor
    print *,message
  end subroutine
'''+includes+'\nend module')
    driver=Path(__file__).with_name('restriction_closure.f90').read_text()
    (out/'closure.f90').write_text(original+driver)


class RestrictionClosure(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_closure(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp); prepare(out)
            for opt in ('-O0','-O2'):
                cmd=['gfortran','-cpp',opt,'-ffp-contract=off','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                     '-Werror=line-truncation','-I'+tmp,str(ROOT/'src/kind.f90'),
                     str(ROOT/'src/parallel_block_scalar_storage.f90'),'closure.f90','-o','check']
                r=subprocess.run(cmd,cwd=out,text=True,capture_output=True)
                self.assertEqual(r.returncode,0,r.stdout+r.stderr)
                for mode in ('valid','poison','nan','extent','full'):
                    r=subprocess.run(['./check',mode],cwd=out,text=True,capture_output=True)
                    if mode=='valid':
                        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
                        self.assertIn('PASS restriction closure',r.stdout)
                    else:
                        self.assertNotEqual(r.returncode,0,mode)
                        self.assertIn({'poison':'poisoned','nan':'non-finite','extent':'extent',
                                       'full':'compact production'}[mode],r.stdout+r.stderr)

if __name__=='__main__':
    unittest.main()
