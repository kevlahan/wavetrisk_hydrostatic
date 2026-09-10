from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from prepare_memory_census import declaration,generate,specification
from run_memory_census import read_census,summarize


FIXTURE='''module fixture
  use iso_fortran_env, only : real64,int64
  type :: box
    real(real64), allocatable :: values(:)
    real(real64), pointer :: alias(:)=>null()
  end type box
  type(box), allocatable :: objects(:)
  integer, allocatable :: empty(:)
contains
end module fixture
'''


class CensusTests(unittest.TestCase):
    def test_dimensions_continuations_and_alias_exclusion(self):
        fields=declaration('real(real64), dimension(:,:), allocatable :: a,b')
        self.assertEqual([(f.name,f.rank) for f in fields],[('a',2),('b',2)])
        types,roots=specification(FIXTURE)
        source,coverage=generate('fixture',roots,types)
        self.assertIn('fixture.objects.values',coverage['allocations'])
        self.assertNotIn('fixture.objects.alias',coverage['allocations'])
        self.assertIn('fixture.objects.alias (pointer target)',coverage['excluded_pointer_targets'])
        self.assertTrue(all(len(line)<=132 for line in source.splitlines()))

    def test_unknown_owning_type_and_declaration_rejected(self):
        with self.assertRaisesRegex(ValueError,'Unresolved'):
            generate('bad',declaration('type(unknown), allocatable :: a(:)'),{})
        with self.assertRaisesRegex(ValueError,'Unsupported owning'):
            declaration('complex, allocatable :: a(:)')

    @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
    def test_actual_allocated_capacity_nested_empty_and_lower_bounds(self):
        types,roots=specification(FIXTURE);routine,_=generate('fixture',roots,types)
        source=FIXTURE.replace('end module fixture',routine+'\nend module fixture')
        source+='''
program test
  use fixture
  integer :: u
  allocate(objects(-1:0),empty(0))
  allocate(objects(-1)%values(-2:2))
  allocate(objects(0)%values(7))
  open(newunit=u,file='result.txt',status='new')
  call census_fixture(u)
  close(u)
  print *, size(objects,kind=int64)*int(storage_size(objects),int64)/8_int64
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'test.f90').write_text(source)
            subprocess.run(['gfortran','-std=f2008','-Wall','-Wextra','-Werror=line-truncation',
                            '-fcheck=all','test.f90','-o','test'],cwd=root,check=True,capture_output=True)
            result=subprocess.run(['./test'],cwd=root,check=True,capture_output=True,text=True)
            rows={line.split()[1]:int(line.split()[0]) for line in (root/'result.txt').read_text().splitlines()}
        self.assertEqual(rows['fixture.objects.values'],12*8)
        self.assertEqual(rows['fixture.empty'],0)
        self.assertEqual(rows['fixture.objects'],int(result.stdout.strip()))

    def test_incomplete_duplicate_or_negative_samples_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'sample'
            for text in ('sample 1 ready\n', 'sample 1 ready\n-1 a.b\n',
                         'sample 1 ready\n1 a.b\n1 a.b\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):read_census(path,{'a.b'})

    def test_aligned_sum_is_not_sum_of_independent_peaks(self):
        manifest={'scope':'test','excluded':'none','coverage':{'a':{'allocations':['a.b']}}}
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'allocation-rank-0.txt').write_text('sample 1 ready\n100 a.b\nsample 2 end\n10 a.b\n')
            (root/'allocation-rank-1.txt').write_text('sample 1 ready\n1 a.b\nsample 2 end\n100 a.b\n')
            result=summarize(root,manifest,2)
            self.assertEqual(result['largest_aligned_phase_sum']['bytes'],110)
            self.assertEqual(sum(r['peak_counted_bytes'] for r in result['ranks']),200)
            (root/'allocation-rank-1.txt').write_text('sample 1 end\n1 a.b\n')
            result=summarize(root,manifest,2)
            self.assertFalse(result['phase_sequences_align'])
            self.assertNotIn('largest_aligned_phase_sum',result)


if __name__=='__main__':unittest.main()
