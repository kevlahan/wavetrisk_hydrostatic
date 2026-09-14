"""Verify new timer nesting and fail-closed split/work analysis."""
from pathlib import Path
import copy
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from analyze_inner_inverse_194 import split_profile,inner_work

ROOT=Path(__file__).resolve().parents[2]


class SplitInverse194(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'),'GNU Fortran required')
    def test_actual_profiler_children_and_disabled_mode(self):
        source='''program check
use parallel_block_profile_mod
implicit none
integer::i,phases(4)=[DP_SCALAR_INVERSE_KERNEL,DP_INNER_VECTOR_KERNEL,DP_SCALAR_INVERSE_KERNEL,DP_OUTER_KERNEL]
real(8)::residual
detail_enabled=.true.
call detail_reset
call detail_enter(DP_STEP)
do i=1,4
 call detail_enter(DP_INVERSE_KERNEL)
 call detail_enter(phases(i))
 call detail_leave(phases(i))
 call detail_leave(DP_INVERSE_KERNEL)
end do
call detail_leave(DP_STEP)
if(detail_calls(DP_INVERSE_KERNEL)/=4.or.detail_calls(DP_SCALAR_INVERSE_KERNEL)/=2.or. &
 detail_calls(DP_INNER_VECTOR_KERNEL)/=1.or.detail_calls(DP_OUTER_KERNEL)/=1) error stop 'split calls'
residual=detail_time(3,DP_INVERSE_KERNEL)-detail_time(1,DP_INVERSE_KERNEL)- &
 detail_time(3,DP_SCALAR_INVERSE_KERNEL)-detail_time(3,DP_INNER_VECTOR_KERNEL)-detail_time(3,DP_OUTER_KERNEL)
if(abs(residual)>1.0d-7) error stop 'parent conservation'
if(abs(sum(detail_time(1,:))-detail_time(3,DP_STEP))>1.0d-7) error stop 'step conservation'
call detail_reset
detail_enabled=.false.
call detail_enter(DP_INVERSE_KERNEL)
call detail_enter(DP_SCALAR_INVERSE_KERNEL)
call detail_leave(DP_SCALAR_INVERSE_KERNEL)
call detail_leave(DP_INVERSE_KERNEL)
if(any(detail_calls/=0).or.any(detail_time/=0)) error stop 'disabled timer active'
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'driver.f90').write_text(source)
            for opt in ('-O0','-O2'):
                subprocess.run(['gfortran',opt,'-fcheck=all','-ffree-line-length-132',str(ROOT/'src/parallel_block_profile.f90'),
                                'driver.f90','-o','driver'],cwd=p,check=True,capture_output=True)
                subprocess.run([str(p/'driver')],cwd=p,check=True,capture_output=True)

    def test_split_rejects_missing_and_unbalanced_regions(self):
        regions=[]
        for i,seconds,calls in [(21,1.,1),(36,6.1,4),(66,3.,2),(67,2.,1)]:
            regions.append(dict(id=i,inclusive_wall_avg=seconds,self_wall_avg=.1 if i==36 else seconds,
                                self_cpu_avg=.01 if i==36 else seconds,global_calls=calls))
        data={'completed':True,'detail_windows':[{'regions':regions,'step_summaries':[
            {'restart':False,'complete':True,'max_self_conservation_error':0.,'self_avg':[1.]*67}]}]}
        with patch('analyze_inner_inverse_194.analyze',return_value=data):
            self.assertEqual(split_profile(Path('unused'))['whole_window_rank_mean']['scalar']['global_calls'],2)
        for failure in ('missing','calls','time','step'):
            bad=copy.deepcopy(data)
            if failure=='missing':bad['detail_windows'][0]['regions'].pop()
            elif failure=='calls':bad['detail_windows'][0]['regions'][2]['global_calls']=3
            elif failure=='time':bad['detail_windows'][0]['regions'][1]['inclusive_wall_avg']=9
            else:bad['detail_windows'][0]['step_summaries'][0]['complete']=False
            with patch('analyze_inner_inverse_194.analyze',return_value=bad),self.assertRaises(ValueError):split_profile(Path('unused'))

    def test_geometry_work_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'run.log'
            prefix='  inner inverse geometry: prepared-cells cell-updates reads reference-reads = '
            p.write_text(prefix+'4 120 128 3840\n'+prefix+'8 240 256 7680\n')
            self.assertAlmostEqual(inner_work(p)['geometry_read_reduction_percent'],96.6666666667)
            p.write_text(prefix+'0 120 3840 3840\n')
            self.assertEqual(inner_work(p)['geometry_read_reduction_percent'],0)
            for values in ('0 120 100 3840','4 120 130 3840','4 120 128 3839','4 0 0 0'):
                p.write_text(prefix+values+'\n')
                with self.assertRaises(ValueError):inner_work(p)
