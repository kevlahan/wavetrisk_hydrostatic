"""Compare the actual narrow boundary update with the historical full-record path."""
from pathlib import Path
import re,shutil,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
class ScalarBoundaryUpdate(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_full_record_equivalence(self):
        source=(ROOT/'src/parallel_block_mpi.f90').read_text()
        update=re.search(r'^    subroutine fill_retained_boundary_record .*?^    end subroutine fill_retained_boundary_record',source,re.M|re.S).group()
        node=re.search(r'^    subroutine fill_boundary_node .*?^    end subroutine fill_boundary_node',source,re.M|re.S).group()
        names=set(re.findall(r'\bBLOCK_[A-Z_]+\b',update+node))
        constants=[]
        for name in sorted(names-{'BLOCK_BOUNDARY_POISON'}):
            match=re.search(r'^  integer, parameter :: '+name+r' =.*$',source,re.M)
            self.assertIsNotNone(match,name)
            constants.append(match.group())
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for name,text in [('new_update.inc',update),('fill_node.inc',node),('constants.inc','\n'.join(constants))]:
                (out/name).write_text(text+'\n')
            shutil.copy2(Path(__file__).with_name('scalar_boundary_reference.inc'),out/'scalar_boundary_reference.inc')
            for opt in ('-O0','-O2'):
                build=subprocess.run(['gfortran',opt,'-fcheck=all','-ffpe-trap=invalid,zero,overflow','-Werror=line-truncation','-I'+tmp,str(ROOT/'src/kind.f90'),str(ROOT/'src/parallel_block_scalar_storage.f90'),str(Path(__file__).with_name('scalar_boundary_update.f90')),'-o','test'],cwd=out,capture_output=True,text=True)
                self.assertEqual(build.returncode,0,build.stdout+build.stderr)
                run=subprocess.run(['./test'],cwd=out,capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                self.assertIn('PASS boundary update',run.stdout)
if __name__=='__main__':unittest.main()
